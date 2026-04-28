clc
clear
clf
close all

OsID=1
OscType={'Rossler' ,'Lorenz','Chen','HR'};
FMtype = char(OscType(OsID))

tlayout = tiledlayout(2,1, "TileSpacing", "tight");

for ii=1:1:1
    [row,col] =ind2sub([3 3],ii);
    
    name=strcat(FMtype,'-Chain-','Case-',num2str(ii),'-Nstp=',num2str(0.25),'.mat')
    
    load(name);
    ax = nexttile
    hold on
    grid on

    spacing = 50;
    marcadores = {'o', 's', '^', 'v', 'd', 'p', 'h', '*', '+', 'x'};
    colores = lines(10);
    handles = zeros(1,10);
    
    
    for i = 1:10
        x_data = KL(:,1);
        y_data = KL(:,i+1);
        indices = 1:spacing:length(x_data);
        
        % Crear un solo objeto combinando línea y marcadores
        handles(i) = plot(x_data, y_data, '-', ...
            'Color', colores(i,:), ...
            'LineWidth', 2, ...
            'Marker', marcadores{i}, ...
            'MarkerSize', 8, ...
            'MarkerFaceColor', colores(i,:), ...
            'MarkerIndices', indices);
    end
    
    legend_labels = {};
    for i = 1:10
        legend_labels{i} = ['$\lambda_{' num2str(i) '}$'];
    end
    legend(handles, legend_labels, 'Interpreter', 'latex', 'Location', 'eastoutside', 'FontSize', 15)
    
    
    ax = gca;
    ax.FontSize = 30;
    title(strcat('$b_',num2str(col),'\rightarrow  c_',num2str(row),'$'),'interpreter','latex','fontsize',35)
    ylabel('$\lambda_{i}$', 'interpreter', 'latex', 'fontsize', 35)
    
    

    Lmax = max(KL(:,2:end), [], 2);
    x_data = KL(:,1);
    
    nexttile
    hold on
    grid on
    

    signos = sign(Lmax);
    cambios = find(diff(signos) ~= 0);
    inicio = 1;
    
    handles = [];
    labels = {};
    
    for k = 1:length(cambios)
        fin = cambios(k);
        if Lmax(inicio) < 0
            h = plot(x_data(inicio:fin), Lmax(inicio:fin), '-g', 'LineWidth', 2);
            if isempty(find(strcmp(labels, '$\lambda_{max}^{\top} < 0$'), 1))
                handles = [handles, h];
                labels{end+1} = '$\lambda_{max}^{\top} < 0$';
            end
        else
            h = plot(x_data(inicio:fin), Lmax(inicio:fin), '-r', 'LineWidth', 2);
            if isempty(find(strcmp(labels, '$\lambda_{max}^{\top} > 0$'), 1))
                handles = [handles, h];
                labels{end+1} = '$\lambda_{max}^{\top} > 0$';
            end
        end
        inicio = fin + 1;
    end
    

    if inicio <= length(x_data)
        if Lmax(inicio) < 0
            h = plot(x_data(inicio:end), Lmax(inicio:end), '-g', 'LineWidth', 2);
            if isempty(find(strcmp(labels, '$\lambda_{max}^{\top} < 0$'), 1))
                handles = [handles, h];
                labels{end+1} = '$\lambda_{max}^{\top} < 0$';
            end
        else
            h = plot(x_data(inicio:end), Lmax(inicio:end), '-r', 'LineWidth', 2);
            if isempty(find(strcmp(labels, '$\lambda_{max}^{\top} > 0$'), 1))
                handles = [handles, h];
                labels{end+1} = '$\lambda_{max}^{\top} > 0$';
            end
        end
    end
    
    lgd=legend(handles, labels, 'Interpreter', 'latex', 'Location', 'eastoutside', 'FontSize', 25,'NumColumns', 1, ...
        'Box', 'on', ...
        'Orientation', 'horizontal');
   
    hold off

    ax = gca;
    ax.FontSize = 35;
    xlabel("$k_d$",'interpreter','latex','fontsize',40)
    axis([-5 5 -0.5 1 ])
    ylabel('$\lambda_{max}^{\top}$', 'interpreter', 'latex', 'fontsize', 40)
    
    
end

exportgraphics(gcf,'ChainStaticDynamicscRosslerb1c1Lyapunov.eps','ContentType','vector')


clc
clear
clf
close all

OsID=1
OscType={'Rossler' ,'Lorenz','Chen','HR'};
FMtype = char(OscType(OsID))

tlayout = tiledlayout(3,3, "TileSpacing", "tight", "Padding", "tight");

for ii = 1:1:9
    [row, col] = ind2sub([3 3], ii);
    
    name = strcat(FMtype, '-Chain-', 'Case-', num2str(ii), '-Nstp=', num2str(0.25), '.mat');
    load(name);

    Lmax = max(KL(:,2:end), [], 2);
    x_data = KL(:,1);
    [x_data Lmax]
    
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

    title(strcat('$b^{(',num2str(col),')}\rightarrow  c^{(',num2str(row),')}$'),'interpreter','latex','FontSize',20);

    ax = gca;
    ax.FontSize = 25;
    axis([-10 10 -1 1])
end


xlabel(tlayout, "$k_d$", 'interpreter', 'latex', 'fontsize', 35)
ylabel(tlayout, '$\lambda_{max}^{\top}$', 'interpreter', 'latex', 'fontsize', 35)

exportgraphics(gcf,'ChainStaticDynamicscRosslerbcLyapunov.eps','ContentType','vector')


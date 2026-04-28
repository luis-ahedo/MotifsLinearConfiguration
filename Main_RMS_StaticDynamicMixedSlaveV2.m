clc
clear
clf
close all


syms x(t) y(t) z(t)
syms x1(t) y1(t) z1(t)
syms x2(t) y2(t) z2(t)
syms x3(t) y3(t) z3(t)
syms a b c K1 K2 t real
syms h1(t) alfa
X=[x(t);y(t);z(t)];
X1=[x1(t);y1(t);z1(t)];
X2=[x2(t);y2(t);z2(t)];
X3=[x3(t);y3(t);z3(t)];
XAll=[X1 ;X2;X3];
Psym=[a b c];

OscType={'Rossler' ,'Lorenz','Chen','HR'};

OsID=1;
[F,Pval] = OscillatorType(OsID)

n= length(X);
N=n*n;

FX1=subs(F,X,X1);
FX2=subs(F,X,X2);
FX3=subs(F,X,X3);

tlayout = tiledlayout(3,3, "TileSpacing", "tight", "Padding", "tight");

[b c]=SISObc(1,n);
ks=1


for ii=1:9

    ii
    [b2 c2]=SISObc(ii,n);

    eqn1=diff(X1,t) ==FX1;
    eqn2=diff(X2,t) ==FX2+ks*b*c*(X1-X2)+ks*b*c*(X3-X2)-b2*h1(t);
    eqn3=diff(X3,t) ==FX3+ks*b*c*(X2-X3);
    eqnh1=diff(h1,t)==-alfa*h1(t)-K1*c2*(X1-X2);

    eqnS=[eqn1;eqn2;eqn3;eqnh1];
    varS=[XAll; h1 ;];
    [eqns,vars] =reduceDifferentialOrder(eqnS,varS);
    [M,F] = massMatrixForm(eqns,vars);

    Par=[Psym K1 alfa];
    options = odeset('Events',@myEventsFcn);
    odeFunction(F,varS,Par,'File','odeMotifsSD');
    tf=10000;
    tspan = [0 tf];
    %
    RMS123=[];
    for k1=-5:1/128:5

        Parval=[Pval k1 1];

        y0 = rand(1,N+1);
        [tk,y,te,ye,ie] = ode45(@(t,y) odeMotifsSD(t,y,Parval), tspan, y0,options);
        % ie
        x1e=y(:,1);
        y1e=y(:,2);
        z1e=y(:,3);

        x2e=y(:,4);
        y2e=y(:,5);
        z2e=y(:,6);

        x3e=y(:,7);
        y3e=y(:,8);
        z3e=y(:,9);


        Erx12=rms(abs(x2e-x1e));
        Ery12=rms(abs(y2e-y1e));
        Erz12=rms(abs(z2e-z1e));

        Erx13=rms(abs(x3e-x1e));
        Ery13=rms(abs(y3e-y1e));
        Erz13=rms(abs(z3e-z1e));

        Erx23=rms(abs(x2e-x3e));
        Ery23=rms(abs(y2e-y3e));
        Erz23=rms(abs(z2e-z3e));

       

        RMS12=mean([rms(Erx12) rms(Ery12) rms(Erz12)]);

        RMS13=mean([rms(Erx13) rms(Ery13) rms(Erz13)]);

        RMS23=mean([rms(Erx23) rms(Ery23) rms(Erz23)]);
        %
        RMS123=vertcat(RMS123,[k1 RMS12 RMS13 RMS23]);


    end

    ax = nexttile;
    plot(RMS123(:,1),RMS123(:,2),'b')
 
    hold on
    grid on
    axis("auto")
    plot(RMS123(:,1),RMS123(:,3),'k')
    plot(RMS123(:,1),RMS123(:,4),'c')
    % axis([-20 20 0 10])
    [rowi,coli] =ind2sub([3 3],ii);
    
    % title(strcat('$b_',num2str(coli),'\rightarrow  c_',num2str(rowi),'$' ...
    %     ,' ',',$b_',num2str(colj),'\rightarrow  c_',num2str(rowj),'$'),'interpreter','latex')
    
    ax.FontSize = 20;

    title(strcat('$b^{(',num2str(coli),')}\rightarrow  c^{(',num2str(rowi),')}$'),'interpreter','latex','FontSize',25);

      
 
end

xlabel(tlayout,'$k_d$','interpreter','latex','fontsize',25)
ylabel(tlayout,'$<\mathbf{e}_{ij}>$','interpreter','latex','fontsize',25)

exportgraphics(gcf,'MixedSlaveStaticDynamicscRossler.eps','ContentType','vector')

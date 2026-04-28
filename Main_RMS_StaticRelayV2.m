clc
clear
clf
close all

syms x(t) y(t) z(t)
syms x1(t) y1(t) z1(t)
syms x2(t) y2(t) z2(t)
syms x3(t) y3(t) z3(t)
syms a b c K1 K2 t real
X=[x(t);y(t);z(t)];
X1=[x1(t);y1(t);z1(t)];
X2=[x2(t);y2(t);z2(t)];
X3=[x3(t);y3(t);z3(t)];
XAll=[X1 ;X2;X3]
Psym=[a b c];

OscType={'Rossler' ,'Lorenz','Chen','HR'};

OsID=1;
[F,Pval] = OscillatorType(OsID)

n= length(X);

FX1=subs(F,X,X1)
FX2=subs(F,X,X2)
FX3=subs(F,X,X3)

N=n*n;

tlayout = tiledlayout(3,3, "TileSpacing", "tight", "Padding", "tight");

[b c]=SISObc(1,n);

ii=1


for jj=1:9

    [jj,jj]
    [b2 c2]=SISObc(jj,n);
    [b3 c3]=SISObc(jj,n);
    [b1 c1]=SISObc(jj,n);

    %% Chain
    eqn1=diff(X1,t) ==FX1+K1*b1*c1*(X2-X1);
    eqn2=diff(X2,t) ==FX2+K1*b2*c2*(X1-X2)+K1*b2*c2*(X3-X2);
    eqn3=diff(X3,t) ==FX3+K1*b3*c3*(X2-X3);


    eqnS=[eqn1;eqn2;eqn3];
    varS=[XAll];
    [eqns,vars] =reduceDifferentialOrder(eqnS,varS);
    [M,FAll] = massMatrixForm(eqns,vars);

    Par=[Psym K1];
    options = odeset('Events',@myEventsFcn);
    odeFunction(FAll,XAll,Par,'File','odeMotifStatic');
    tf=10000;
    tspan = [0 tf];
    RMS123=[];



    for k1=-5:1/128:5

        Parval=[Pval k1];

        y0 = rand(1,N);
        [tk,y,te,ye,ie] = ode45(@(t,y) odeMotifStatic(t,y,Parval), tspan, y0,options);
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


        Erx12=abs(x2e-x1e);
        Ery12=abs(y2e-y1e);
        Erz12=abs(z2e-z1e);

        Erx13=abs(x3e-x1e);
        Ery13=abs(y3e-y1e);
        Erz13=abs(z3e-z1e);

        Erx23=abs(x2e-x3e);
        Ery23=abs(y2e-y3e);
        Erz23=abs(z2e-z3e);

        RMS12=mean([rms(Erx12) rms(Ery12) rms(Erz12)]);

        RMS13=mean([rms(Erx13) rms(Ery13) rms(Erz13)]);

        RMS23=mean([rms(Erx23) rms(Ery23) rms(Erz23)]);

        RMS123=vertcat(RMS123,[k1 RMS12 RMS13 RMS23]);

    end

    ax = nexttile 
    plot(RMS123(:,1),RMS123(:,2),'b')
    hold on
    grid on
     ax.FontSize = 20;

    plot(RMS123(:,1),RMS123(:,3),'k')
    plot(RMS123(:,1),RMS123(:,4),'c')

    [rowi,coli] =ind2sub([3 3],jj);
    [rowj,colj] =ind2sub([3 3],jj);

    % title(strcat('$b_',num2str(2),'^',num2str(coli),'\rightarrow  c_',num2str(2),'^',num2str(rowi),'$' ...
        % ,' ',',$b_',num2str(3),'^',num2str(colj),'\rightarrow  c_',num2str(3),'^',num2str(rowj),'$'),'interpreter','latex','FontSize',20)

title(strcat('$b_k','^{(',num2str(colj),')}\rightarrow  c_k','^{(',num2str(rowj),')}$'),'interpreter','latex','FontSize',25)
        
   

end

xlabel(tlayout,'$k$','interpreter','latex','fontsize',25)
ylabel(tlayout,'$<\mathbf{e}_{ij}>$','interpreter','latex','fontsize',25)

exportgraphics(gcf,'RelayStaticRossler.eps','ContentType','vector')

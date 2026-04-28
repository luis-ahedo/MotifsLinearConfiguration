clc
clear
clf
close all

syms x(t) y(t) z(t)
syms x1(t) y1(t) z1(t)
syms x2(t) y2(t) z2(t)
syms a b c kd ks t real
syms h(t) alfa

X=[x(t);y(t);z(t)];
X1=[x1(t);y1(t);z1(t)];
X2=[x2(t);y2(t);z2(t)];

n= length(X);
N=n*n+1;
syms Y(t) [N 1]

H=[h];
Psym=[a b c];

OscType={'Rossler' ,'Lorenz','Chen','HR'};
OsID=1;
FMtype=char(OscType(OsID));
[F,Pval] = OscillatorType(OsID);
F=subs(F,Psym,Pval);


[bs cs]=SISObc(1,n);
ks=2;

for ii=1:1

    [b2 c2]=SISObc(ii,n);

    eqn=diff(X,t)==F;
    eqnFx=[eqn];
    varFx=[X ];

    eqnDN=F;
    eqnHDN=-alfa*h(t) -kd*c2*(X1-X2);

    [eqns,vars] =reduceDifferentialOrder(eqnHDN,H);
    DFH=jacobian(eqns,vars);

    [eqns,vars] =reduceDifferentialOrder(eqnDN,X);
    DFX=jacobian(eqns,vars);

    M1=[DFX-ks*bs*cs zeros(n,n) zeros(n,n) b2];
    M2=[ks*bs*cs DFX-ks*bs*cs zeros(n,n) -b2];
    M3=[ks*bs*cs zeros(n,n) DFX-ks*bs*cs zeros(n,1)];
    M4=[-kd*c2 zeros(1,n) zeros(1,n) DFH ];

    eqnY=diff(Y,t)==vertcat(M1,M2,M3,M4)*Y;

    eqnFxY=[eqnFx; eqnY];
    varFxY= [varFx; Y];

    [eqns,vars] =reduceDifferentialOrder(eqnFxY,varFxY);
    [M,FXY] = massMatrixForm(eqns,vars);

    Par=[alfa kd];
    odeFunction(FXY,varFxY,Par,'File','odeVariational');

    KD=-10:1/64:10;
    m=length(KD);

    KL=[];

    % clc
    % tic
    % Parval=[1  KD];
    % LyaExTrans= MaxLyaExp(Parval,vars)
    % toc

    for i=1:m
        Parval=[1  KD(i)]
        LyaExTrans= MaxLyaExp(Parval,N);
        KL=vertcat(KL,[KD(i) LyaExTrans])

    end

    name=strcat(FMtype,'-Chain-','Case-',num2str(ii),'-Nstp=',num2str(0.25),'.mat')
    save(name,'KL')

end


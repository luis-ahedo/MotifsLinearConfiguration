
function [FM,OscPar,FMtype,Range] = OscillatorType(idOsc)

OscType={'Rossler' ,'Lorenz','Chen','HR'};
FMtype = char(OscType(idOsc));

syms a b c x(t) y(t) z(t)
Xm=[x(t);y(t);z(t)];

switch FMtype
    case 'Rossler'
        FM=[ -y(t)-z(t); x(t)+a*y(t); b+z(t)*(x(t)-c)];
        ap=0.2;
        bp=0.2;
        cp=5.7;
    case 'Lorenz'
        FM=[ a*(y(t)-x(t)); x(t)*(b-z(t))-y(t); x(t)*y(t)-c*z(t)];
        ap=10;
        bp=28;
        cp=2;
    case 'Chen'
        FM=[ a*(y(t)-x(t));(b-a-z(t))*x(t)+b*y(t);x(t)*y(t)-c*z(t)]
        ap=35;
        bp=28;
        cp=8/3;
    case 'HR'
        FM=[ ym(t)+3*(xm(t))^2-(xm(t))^3-zm(t)+a;1-5*(xm(t))^2-ym(t);-b*zm(t)+b*c*(xm(t)+1.6)]
        ap=3.2;
        bp=0.006;
        cp=4;
    otherwise
        FM=[ -ym(t)-zm(t); xm(t)+a*ym(t); b+zm(t)*(xm(t)-c)]
end

OscPar=[ap,bp,cp];
end
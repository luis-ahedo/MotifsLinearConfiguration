function LyaExTrans= MaxLyaExp(Parval,N)
N=N;
tf=10000;
tspan=[0 tf];
dt=0.001;
t  = (tspan(1):dt:tspan(end)).';
n=3;
x0=rand(3,1);
X0=repmat(x0, 1, N);

LyaEx= zeros(1,N);
Y0=rand(N,N);           %%% IC variational
[LyaEx, Y0]=GramSchmidtV1(LyaEx,Y0);

x0=vertcat(X0,Y0);

ns = length(t);         %%% Number of samples
xk = x0;                %%% Current state
hdt = 0.5 * dt;         %%% mean step

for k = 1:ns
    % ns-1
    % k
    % From k to k+1
    tk = t(k);                                    % Current time
    k1 = odeVariational(tk,xk,Parval);                           % RK derivatives
    k2 = odeVariational(tk + hdt, xk + hdt * k1,Parval);                % ...
    k3 = odeVariational(tk + hdt, xk + hdt * k2,Parval);
    k4 = odeVariational(tk +  dt, xk +  dt * k3,Parval);
    xk = xk + (k1 + 2.*k2 + 2.*k3 + k4).*(dt./6); % Updated state
    
    Y0=xk(n+1:end,1:N);
    [LyaEx, Y0]=GramSchmidtV1(LyaEx,Y0);
    xk(n+1:end,1:N)=Y0;
    mask = isnan(LyaEx) ;
    hasAny = any(mask(:));
    if (hasAny==1)
        tk
        LyaEx
        break;
    end
    
end

LyaExTrans=(LyaEx/tf);

end
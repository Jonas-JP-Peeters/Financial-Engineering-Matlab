%% Calculate the price of a European Option (Simple Greek)
T = 3; S = 100; %Time to Maturity (T) and Stock Price (S)
r = 0.0322; sigma = 0.1815; %RF-Rate (r) and volatility (sigma)
X = 100; q = 0.033; %Strike Price (X) and Dividend Yield (q)
d1 = (log(S/X)+(r-q+sigma^2/2)*T)/(sigma*sqrt(T));
d2 = d1-sigma*sqrt(T);
C = (S*exp((r-q)*T)*normcdf(d1)-X*normcdf(d2))*exp(-r*T);
P = (X*normcdf(-d2)-S*exp((r-q)*T)*normcdf(-d1))*exp(-r*T);
%% Calculate the price of a European Option (Binomial)
T = 6/12; dt = 3/12; N = T/dt+1;
X = 51; S0 = 50; S = triu(ones(N,N))*S0;
C = triu(ones(N,N)); P = triu(ones(N,N));
sigma = 0.4; r = 0.05; q = 0.00;
%u = exp(sigma*sqrt(dt)); d = 1/u; 
u = 1.06; d = 0.95;
p = (exp((r-q)*dt)-d)/(u-d);
for i = N:-1:1
    for j = 1:i
        S(j,i) = S0*u^(i-j)*d^(j-1);
        if i == N
            C(j,i) = max(S(j,i)-X,0);
            P(j,i) = max(X-S(j,i),0);
        else
            C(j,i) = (C(j,i+1)*p+C(j+1,i+1)*(1-p))*exp(-r*dt);
            P(j,i) = (P(j,i+1)*p+P(j+1,i+1)*(1-p))*exp(-r*dt);
        end
    end
end
%% Calculate the price of a European Option (Cox, Ross, Rubinstein)
dt = 3/12; T = 6/12; r = 0.05;
n = T/dt; %Number of discrete steps (n)
%u = sigma*sqrt(T/n); d = 1/u;
u = 1.06; d = 0.95;
S = 50; X = 51; syms k;
r = (r+1)^(dt); p = (r-d)/(u-d);
C = double(...
        symsum(...
            nchoosek(n,k)*p^k*(1-p)^(n-k)*...
            piecewise(u^k*d^(n-k)*S-X > 0, (u^k*d^(n-k)*S-X), 0),...
            k, 0, n)...
        /r^n);
P = double(...
        symsum(...
            nchoosek(n,k)*p^k*(1-p)^(n-k)*...
            piecewise(X-u^k*d^(n-k)*S > 0, (X-u^k*d^(n-k)*S), 0),...
            k, 0, n)...
        /r^n);
%% Calculate the price of a European Option (Simulation)
Q = 1000; % Number of simulations (Q)
X = 300; S0 = 300; S = ones(1,Q)*S0;
sigma = 0.2; T = 2/12; r = 0.08; q = 0.03;
E = rand([1 Q]); E = norminv(E); % Errors
S(:,:) = S(:,:).*exp((r-q-0.5*sigma^2)*T+sigma.*sqrt(T)*E(:,:));
C = max(S(:,:)-X,0); C = mean(C)*exp(-r*T);
P = max(X-S(:,:),0); P = mean(P)*exp(-r*T);
%% Calculate the price of a American Option (Binomial)
N = 5; % Number of steps (N)
X = 51; S0 = 50; S = triu(ones(N,N))*S0;
C = triu(ones(N,N)); P = triu(ones(N,N));
sigma = 0.4; dt = 3/12; r = 0.05; q = 0.00; T = 6/12;
u = exp(sigma*sqrt(dt)); d = 1/u; p = (exp((r-q)*dt)-d)/(u-d);
for i = N:-1:1
    for j = 1:i
        S(j,i) = S0*u^(i-j)*d^(j-1);
        if i == N
            C(j,i) = max(S(j,i)-X,0);
            P(j,i) = max(X-S(j,i),0);
        else
            C(j,i) = max(S(j,i)-X,...
                (C(j,i+1)*p+C(j+1,i+1)*(1-p))*exp(-r*dt));
            P(j,i) = max(X-S(j,i),...
                (P(j,i+1)*p+P(j+1,i+1)*(1-p))*exp(-r*dt));
        end
    end
end
%% Calculate the price of a Asian Option (Simulation)
N = 36; Q = 100; % Number of steps (N) % Number of simulations (Q)
S0 = 1; S = ones(1,Q,N)*S0;
E = rand([1 Q N]); E = norminv(E); % Errors
sigma = 0.2; dt = 1/12; r = 0.1; q = 0.00;
for i = 2:N
    S(:,:,i) = S(:,:,i-1).*exp((r-q-0.5*sigma.^2)*dt+sigma.*(sqrt(dt)*E(:,:,i)));
end
A = mean(mean(S,3))*exp(-r*T);
%% Calculate the price of a Quanto Option (Simple Greek)
T = 3; S = 100; %Time to Maturity (T) and Stock Price (S)
r_drift = 0.0322; sigma = 0.1815; %RF-Rate (r) and volatility (sigma)
X = 100; q = 0.033; %Strike Price (X) and Dividend Yield (q)
d1 = (log(S/X)+(r_drift-q-c+sigma^2/2)*T)/(sigma*sqrt(T));
d2 = d1-sigma*sqrt(T); 
corr = -0.3; vol_fx = 0.09;
r_payoff = 0.0442;
c = corr*sigma*vol_fx;
C = (S*exp((r_drift-q-c)*T)*normcdf(d1)-X*normcdf(d2))*exp(-r_payoff*T);
P = (X*normcdf(-d2)-S*exp((r_drift-q-c)*T)*normcdf(-d1))*exp(-r_payoff*T);
%% Calculate the various Greeks of an Option
S = 200; X = 100;
r = 0.04; q = 0.00; sigma = 0.2; T = 5;
% Delta
[CallDelta,PutDelta] = blsdelta(S,X,r,T,sigma,q);
% Gamma
Gamma = blsgamma(S,X,r,T,sigma,q);
% Vega
Vega = blsvega(S,X,r,T,sigma,q);
% Rho
[CallRho,PutRho] = blsrho(S,X,r,T,sigma,q);
% Theta
[CallTheta,PutTheta] = blstheta(S,X,r,T,sigma,q);
%% Delta, Gamma and Vega Hedge a Portfolio
O = [-1000 0.5 2.2 1.8;
    -500 0.8 0.6 0.2;
    -2000 -0.4 1.3 0.7;
    -500 0.7 1.8 1.4];
O = O(:,2:4).*O(:,1);
O = sum(O);
H = [1 0 0;
    0.6 1.5 0.8;
    0.1 0.5 0.6];
syms x y z
params = [x y z];
eqs = H'*params' == -O';
solx = solve(eqs,[x y z]);
display([solx.x solx.y solx.z])
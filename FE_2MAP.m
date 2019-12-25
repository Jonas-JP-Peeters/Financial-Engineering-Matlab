%% Model a 'normal' Brownian Motion
N = 100; % Number of steps
X = zeros(1,N);
E = rand([1 N]); E = norminv(E); % Errors chosen from normal distribution
a = 0.05; b = 0.1; dt = 1/12; % Parameters (a = drift, b = vol)
for i = 2:N
    X(i) = X(i-1)+ a*dt + b*sqrt(dt)*E(i);
end
%% Model a 'geometric' Brownian Motion
X2 = ones(1,N);
% Use same errors and parameters as the 'normal' Brownian motion
for i = 2:N
    X2(i) = X2(i-1)+ X2(i-1)*(a*dt + b*sqrt(dt)*E(i));
end
%% Model a mean-reverting Brownian Motion
X3 = zeros(1,N);
theta = 0.1; mean = 0;
% Use same errors and parameters as the 'normal' Brownian motion
for i = 2:N
    X3(i) = X3(i-1) + theta*(mean - X3(i-1))*dt + b*sqrt(dt)*E(i);
end
%% Model a Realistic Stock
X4 = ones(1,N);
mean = 0.05;
for i = 2:N
    X4(i) = X4(i-1)*exp((mean-0.5*b^2)*dt+b*sqrt(dt)*E(i));
end
%% Calculate the probability distribution of the stock price St
S0 = 50; r = 0.10; sigma = 0.3; dt = 1/252; 
% Calculate the expected stock price on time t
eSt = S0*exp((r-sigma^2/2)*dt);
% Calculate the confidence interval
cl = 0.66; lb = norminv((1-cl)/2); ub = -lb;
CI = [S0*exp((r-sigma^2/2)*dt+lb*sigma*sqrt(dt))...
    S0*exp((r-sigma^2/2)*dt+ub*sigma*sqrt(dt))];
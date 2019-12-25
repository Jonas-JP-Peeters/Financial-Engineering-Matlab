%% Model a correlated geometric Brownian Motion
N = 100; % Number of steps
X = ones(2,N); % Two rows bvs two stocks
E = rand([2 N]); E = norminv(E); % Errors
a = [0.05;0.10]; b = [0.15;0.30]; dt = 1/12; % Parameters
C = eye(2); C(C==0) = 0.80; C = chol(C,'lower');
for i = 2:N
    X(:,i) = X(:,i-1).*(1 + a*dt + b.*(sqrt(dt)*C*E(:,i)));
end
%% Simulation of correlated stock prices
N = 100; % Number of steps
Q = 100; % Number of simulations
X = ones(2,Q,N); % Two rows bvs two stocks
E = rand([2 Q N]); E = norminv(E); % Errors
mean = [0.05;0.10]; b = [0.15;0.30]; dt = 1/12; % Parameters
C = eye(2); C(C==0) = 0.30; C = chol(C,'lower');
for i = 2:N
    X(:,:,i) = X(:,:,i-1).*exp((mean-0.5*b.^2)*dt+b.*(sqrt(dt)*C*E(:,:,i)));
end
plot(permute(X(:,1,:),[3,1,2])) % Plots the first simulation
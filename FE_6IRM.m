%% Calculate the Macaulay's Duration of a Bond (Yield to Maturity)
T = 20; c = 0.03; y = 0.5; % Coupon Rate (c) and # Years (T)
BOND = [(1:T)' ones(T,1) (1:T)' ones(T,1) ones(T,1)];
BOND(1:(T-1),2) = BOND(1:(T-1),2)*c;
BOND(T,2) = BOND(T,2)+c;
BOND(:,3) = BOND(:,2)./(1+y).^BOND(:,1); % PV of Payments
BOND(:,4) = BOND(:,3)/sum(BOND(:,3)); % Weights of Payments
BOND(:,5) = BOND(:,1).*BOND(:,4); % Weighted Payments
D = sum(BOND(:,5));
%% Calculate Forward Rates from Spot Rates
S = [0.0311 0.0333 0.0344 0.0353 0.0359...
    0.0364 0.0369 0.0374 0.0378 0.0382];
F = zeros(size(S)); S = S+1;
for i = 1:size(S,2)
    F(i) = S(i)^i/prod(F(1:i-1));
end
F = F-1;
%% Fit Forward Rates using a Functional Form
delta1 = 0.015; delta2 = -0.002; delta3 = -0.01; lambda = 0.9;
t = 0; T = 1:10; F = ones(1,10);
F = delta1+delta2*exp(-lambda*(T-t))+delta3*lambda*exp(-lambda*(T-t));
%% Calculate Forward Rate Volatility from Spot Rate Volatility
s = [0.1136 0.1537 0.1678 0.1730 0.1740...
    0.1735 0.1722 0.1706 0.1690 0.1672];
sigma = zeros(size(s));
for i = 1:size(s,2)
    sigma(i) = sqrt(i*s(i)^2-sum(sigma(1:i-1).^2));
end
%% Fit Forward Rate Volatility using a Functional Form
a = 0.01; b = 0.1; c = 0.2; d = 0.05;
t = 0; T = 1:10; sigma = ones(1,size(T,2));
Sigma = (a+b*(T-t)).*exp(-c*(T-t))+d;
%% Fit Forward Rate Correlation using a Functional Form
alpha = 0.01; beta = 0.15;
t = 0; T = 1:10 ; rho = ones(size(T,2),size(T,2));
for i = 1:size(T,2)
    rho(:,i) = alpha + (1-alpha)*exp(-beta*abs(T-T(i)));
end
%% LIBOR Market Model
sigma = Sigma; Q = 5; C = chol(rho,'lower');
LMM = zeros(size(F,2),Q,size(F,2)); LMM(:,:,1) = repmat(F',1,Q); 
R = zeros(size(F,2),Q); VVC = sigma'*sigma.*rho;
DT1 = zeros(size(F,2),Q); DT2 = zeros(size(F,2),Q);
for i = 2:size(F,2)
    q = size(F,2)-i+1; sigma = sigma(1:q);
    VVC = VVC(1:q,1:q); C = C(1:q,1:q);
    E = rand([q Q]); E = norminv(E);
    VT = sigma'.*(C*E); R = R(1:q,1:Q);
    DT1 = DT1(1:q,1:Q); DT2 = DT2(1:q,1:Q); 
    for j = 1:q
        k_temp = zeros(j,Q);
        for k = 1:j
            k_temp(k,:) = VVC(k,j)*(LMM(k+i-1,:,i-1)/(1+LMM(k+i-1,:,i-1)));
        end
        DT1(j,:) = sum(k_temp);
        R(j,:) = LMM(j+i-1,:,i-1).*exp((DT1(j,:)-0.5*sigma(j).^2)+VT(j,:));
        for k = 1:j
            k_temp(k,:) = VVC(k,j)*(R(k,:)/(1+R(k,:)));
        end
        DT2(j,:) = sum(k_temp);
        LMM(j+i-1,:,i) = LMM(j+i-1,:,i-1).*exp((0.5*(DT1(j,:)+DT2(j,:))-0.5*sigma(j).^2)+VT(j,:));
    end  
end
p = 0; YC = zeros(size(F,2)-p,Q); % Yield curve after 'p' years
for i = 1+p:size(F,2)
    YC(i-p,:) = (prod(LMM(1+p:i,:,1+p)+1,1)).^(1/i);
end
plot(YC)
%% Calculate the Option Package Value of a Product
%Each row of the matrix O contains an option
%O = [1/0 (Call or Put),
%     1/-1 (Buy or Sell),
%     0-1 (Participation),
%     'double' (Strike Price)];
O = [1, -1, 1, 160;
    0, -1, 1, 100;
    0, 1, 1, 90];
%Parameters
S = 100; r = 0.035; q = 0.02; sigma = 0.2; T = 4;
BL = ones(size(O,1),2); OP = ones(size(O,1),1);
for i = 1:size(O,1)
    [C,P] = blsprice(S,O(i,4),r,T,sigma,q);
    BL(i,1) = C; BL(i,2) = P;
    if O(i,1) == 1
        OP(i,1) = BL(i,1)*O(i,2)*O(i,3);
    else
        OP(i,1) = BL(i,2)*O(i,2)*O(i,3);
    end
end
OPS = sum(OP);
%% Valuate the Price of a Timing Optimizer
X = 100; r = 0.035; q = 0.02; sigma = 0.2; T = 4;
N = T*12; Q = 1000; m = 6;
S0 = 100; S = ones(Q,N)*S0; dt = 1/12;
E = rand([Q N]); E = norminv(E); % Errors
for i = 2:N
    S(:,i) = S(:,i-1).*exp((r-q-0.5*sigma.^2)*dt+sigma.*(sqrt(dt)*E(:,i)));
end
OPT = max(S(:,N)-min(X,min(S(:,2:(m+1)),[],2)),0);
TO = mean(OPT)*exp(-r*T);
%% Valuate the Price of an Autocallable
S = 100; X = 90; r = 0.035; q = 0.02; sigma = 0.2; T = 3;
c = 0.1; dt = 1; N = T/dt+1; Q = 10000; % Coupon Rate (c)
S0 = 100; S = ones(Q,N)*S0; 
PYO = zeros(Q,N-1); PYS = zeros(1,N-1); STOP = zeros(Q,1);
E = rand([Q N]); E = norminv(E); % Errors
for i = 2:N
    S(:,i) = S(:,i-1).*exp((r-q-0.5*sigma.^2)*dt+sigma.*(sqrt(dt)*E(:,i)));
    penalty = (i == N & S(:,i) < X & STOP ~= 1);
    coupon = (S(:,i) > X & STOP ~= 1);
    PYO(:,i-1) = penalty.*(-(X-S(:,i))/S0) + coupon.*(c*(i-1));
    STOP = STOP+coupon;
    PYS(1,i-1) = mean(PYO(:,i-1))*exp(-r*(i-1));
end
AC = sum(PYS);
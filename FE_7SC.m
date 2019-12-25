%% Calculate the Implied Probability of Default using the Premium
r = 0.01; p = 0.04; T_m = 5;
PD = 1 - ((1+r)/(1+r+p))^T_m;
%% Calculate the Conditional Probability of X Defaults.
S = [0.002384 0.001225 0.005467 0.004275 0.0035 0.003933];
X = [0.25 0.25 0.25 0.25 0.25 0.25]; t = 0; T_m = 1:5;
% Unconditional Probability of Default
UPD = 1-exp(-T_m'*(S./(1-X)));
dm = 1/2500; m = dm:dm:(1-dm); m = norminv(m);
% Conditional Probability of Default
T = 5; c2 = 0.7;
CPD = normcdf((norminv(UPD(T,:))-m'*sqrt(c2))/sqrt(1-c2));
DF_0 = prod(1-CPD,2); DF_0 = sum((DF_0(1:end-1)+DF_0(2:end))*dm/2);
DF = zeros(size(m,2),size(S,2),size(S,2)+1);
DF(:,:,1) = CPD; q = size(S,2); DF_s = [DF_0 zeros(1,q)];
for i = q:-1:1
    for j = (q-i+1):q
        if j == (q-i+1)
            DF(:,j,q+2-i) = prod(DF(:,1:(q+1-i),1),2);
        else
            if i == q
                DF(:,j,q+2-i) = DF(:,j-1,q+2-i).*(1-DF(:,j,1))+...
                    DF(:,j,1).*(1-DF(:,j-1,q+2-i));
            else
                DF(:,j,q+2-i) = DF(:,j-1,q+2-i).*(1-DF(:,j,1))+...
                    DF(:,j,1).*DF(:,j-1,q+1-i);
            end
        end
    end
    DF_s(q+2-i) = sum((DF(1:end-1,q,q+2-i)+DF(2:end,q,q+2-i))*dm/2);
end
%% Calculate the Coupon of a First-to-Default
load('DF_S'); r = 0.01; c2_bm = 0.7; T = 1:size(DF_S,1); X = 0.25;
FV5 = DF_S(size(DF_S,1),1,int8(c2_bm*10))*exp(-size(DF_S,1)*r);
CEV = [1 DF_S(:,1,int8(c2_bm*10))']-[DF_S(:,1,int8(c2_bm*10))' 0];
CEV = sum(CEV(1:5).*exp(-T*r)*(1-X));
c = (1-FV5+CEV)/sum(DF_S(:,1,int8(c2_bm*10)).*exp(-T'*r));
P = zeros(size(DF_S,3),1)';
for i = 0.1:0.1:1
    FV5 = DF_S(size(DF_S,1),1,int8(i*10))*exp(-size(DF_S,1)*r);
    CEV = [1 DF_S(:,1,int8(i*10))']-[DF_S(:,1,int8(i*10))' 0];
    CEV = sum(CEV(1:5).*exp(-T*r)*(1-X));
    P(int8(i*10)) = FV5 - CEV + c*sum(DF_S(:,1,int8(i*10)).*exp(-T'*r));
end
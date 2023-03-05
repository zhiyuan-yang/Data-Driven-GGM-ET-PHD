function [C_k1] = predictCPHDcardinalityDistribution(C_k,model)

max_card = model.max_card;
card_gam_k = model.birth.card;
p_S_k = model.kin.pS;

% Allocate memory
C_k1 = zeros(size(C_k));

% Precompute binomial coefficients
BinCoeff = binomial_coeff_matrix(0:max_card,0:max_card);

sum_over_l = zeros(1,max_card+1);
for j = 0:max_card
    L = j:max_card;
    LJ = L-j;
    sum_over_l(j+1) = sum(BinCoeff(L+1,j+1)'.*C_k(L+1).*(p_S_k^j).*((1-p_S_k).^LJ),2); % L+1,j+1¿¼ÂÇ¾ØÕó´Ó1¿ªÊ¼
    C_k1(j+1) = sum(card_gam_k(j-(0:j)+1).*sum_over_l((0:j)+1));
end

function [BC] = binomial_coeff_matrix(N,K)

% Computes the binomial coefficient n!/((n-k)!k!) for the combinations of n
% in N and k in K

n_N = length(N);
n_K = length(K);

NN = repmat(N(:),1,n_K);
KK = repmat(K(:)',n_N,1);
NNKK = max(NN-KK,0);

logBC = log(factorial(NN))-log(factorial(NNKK))-log(factorial(KK));
BC = exp(logBC);
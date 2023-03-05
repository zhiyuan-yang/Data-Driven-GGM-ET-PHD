function [xPts,w_m,w_c] = UT(x,P,alpha,beta,kappa)
% UT transform
n = size(x,1);
lamda = alpha^2 * (n + kappa) - n;

% Calculate Weights Mean
w_m(1) = lamda / (n + lamda);
w_m(2 : 2 * n + 1) = 1 / (2 * (n + lamda));
% Calculate Weights Covariance
w_c(1) = (lamda / (n + lamda)) + (1 - alpha^2 + beta);
w_c(2 : 2 * n + 1) = 1 / (2 * (n + lamda));

%Calculate Sigma Points
A = sqrt(n + lamda) * chol(P)';
xSigma = [zeros(size(x)) A -A];
xPts = xSigma + repmat(x, 1, size(xSigma, 2));

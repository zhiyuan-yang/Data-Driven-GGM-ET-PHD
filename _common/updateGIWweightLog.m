function [Lik] = updateGIWweightLog(gam,p_D,Beta_FA,d,Nm,S,Xpred,Xcorr,nupred,nucorr)

% Calculates the constant for weight update in the Gaussian Inverse Wishart
% Random Matrix PHD filter.
%
% gam - Expected number of measurements
% p_D - probability of detection
% Beta_FA - Clutter rate, (per scan per volume)
% d - Dimension of extension
% Nm - Number of measurements
% S - Scalar innovation factor
% Xpred - Predicted estimate of extension
% Xcorr - Corrected estimate of extension
% nupred - Predicted degrees of freedom for Wishart distribution
% nucorr - Corrected degrees of freedom for Wishart distribution

logC = -gam+log(p_D)-(d/2)*log(Nm*S);
logD = (nupred/2)*log(det(Xpred))-(nucorr/2)*log(det(Xcorr));
logP = Nm*(log(gam)-log(Beta_FA)-(d/2)*log(pi));
logG = gammaln(nucorr/2)+gammaln((nucorr-1)/2)-gammaln(nupred/2)-gammaln((nupred-1)/2);
logL = logC+logD+logP+logG;
Lik = exp(min(logL,709)); %To avoid Lik=inf
end
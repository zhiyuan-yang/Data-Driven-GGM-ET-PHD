function [Lik] = updateGammaweightLog(a_kW,a_kk,b_kW,b_kk)
logL = gammaln(a_kW) - gammaln(a_kk) + a_kk*log(b_kk) - a_kW*log(b_kW);
Lik = exp(logL);

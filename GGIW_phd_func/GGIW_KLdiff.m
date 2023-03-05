function [KLdiff,w_diff,G_diff,N_diff,IW_diff] = ...
    GGIW_KLdiff(w1,a1,b1,m1,P1,v1,V1,w2,a2,b2,m2,P2,v2,V2)

[KLdiv_w1,KLdiv_G1,KLdiv_N1,KLdiv_IW1] = GGIW_KLdiv(w1,a1,b1,m1,P1,v1,V1,w2,a2,b2,m2,P2,v2,V2);
[KLdiv_w2,KLdiv_G2,KLdiv_N2,KLdiv_IW2] = GGIW_KLdiv(w2,a2,b2,m2,P2,v2,V2,w1,a1,b1,m1,P1,v1,V1);

w_diff = KLdiv_w1+KLdiv_w2;
G_diff = KLdiv_G1+KLdiv_G2;
N_diff = KLdiv_N1+KLdiv_N2;
IW_diff = KLdiv_IW1+KLdiv_IW2;

KLdiff = w_diff+G_diff+N_diff+IW_diff;
end

function [KLdiv_w,KLdiv_G,KLdiv_N,KLdiv_IW] = GGIW_KLdiv(w1,a1,b1,m1,P1,v1,V1,w2,a2,b2,m2,P2,v2,V2)

% Function that computes the KL-divergence between two
% Gamma-Gaussian-inverse Wishart (GGIW) distributions
%
% p(g,x,X) = Gam(g;a1,b1)N(x;m1,P1)IW(X;v1,V1)
% q(g,x,X) = Gam(g;a2,b2)N(x;m2,P2)IW(X;v2,V2)
%
% D_KL = D(p||q) = int p ln(p/q) dx

% Dimension of extension
d = size(V1,1);
% Dimension of kinematical state
n_x = length(m1);

KLdiv_w = w1*(log(w1)-log(w2));

KLdiv_G = a1*log(b1)-a2*log(b2)+gammaln(a2)-gammaln(a1)...
    +(a1-a2)*(psi(0,a1)-log(b1))+a1*(b2/b1-1);

KLdiv_N = w1*(...
    -0.5*log(det(P1))+0.5*log(det(P2))...
    -0.5*n_x+0.5*(m1-m2)'*(P2\(m1-m2))...
    +0.5*trace(P2\P1)...
    );

KLdiv_IW = w1*(...
    0.5*(v1-d-1)*log(det(V1))-0.5*(v2-d-1)*log(det(V2))...
    +sum(gammaln((v2-d-(1:d))/2)-gammaln((v1-d-(1:d))/2))...
    +0.5*(v2-v1)*(log(det(V1))-sum(psi(0,(v1-d-(1:d))/2)))...
    +trace(-0.5*(v1-d-1)*(V1\(V1-V2)))...
    );
end
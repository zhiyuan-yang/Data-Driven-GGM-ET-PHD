function [w_k,a_k,b_k,m_k,P_k,nu_k,V_k,J_k] = ...
    GGIW_PHD_update(w_kk1,a_kk1,b_kk1,m_kk1,P_kk1,nu_kk1,V_kk1,J_kk1,Zk,model)

H_k = model.GGIW.H_k;
% Probabilities of survival and detection
p_D_k = model.obs.pD;

% Identity matrix same size as extension
d = model.GGIW.d;
Id = model.GGIW.Id;
% Clutter rate
beta_FA = model.obs.lambda_c * model.obs.pdf_c;

% Allocate some memory

nZk = size(Zk,2);
w_k = zeros(1,(nZk+1)*J_kk1);
a_k = zeros(1,(nZk+1)*J_kk1);
b_k = zeros(1,(nZk+1)*J_kk1);
m_k = zeros(size(m_kk1,1),(nZk+1)*J_kk1);
P_k = zeros(size(P_kk1,1),size(P_kk1,2),(nZk+1)*J_kk1);
nu_k = zeros(1,(nZk+1)*J_kk1);
V_k = zeros(size(V_kk1,1),size(V_kk1,2),(nZk+1)*J_kk1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: The expected number of generated measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gam(1:J_kk1) = a_kk1(1:J_kk1)./b_kk1(1:J_kk1);   %beta_D;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: Update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate memory for probability of detection
pD = p_D_k*ones(1,J_kk1);
for j = 1:J_kk1
    % Gaussian component corresponding to no detections
    w_k(1,j) = (1-(1-exp(-gam(j)))*pD(j))*w_kk1(1,j);
    a_k(1,j) = a_kk1(1,j);
    b_k(1,j) = b_kk1(1,j);
    m_k(:,j) = m_kk1(:,j);
    P_k(:,:,j) = P_kk1(:,:,j);
    nu_k(1,j) = nu_kk1(1,j);
    V_k(:,:,j) = V_kk1(:,:,j);
end

% l is an index that keeps track of how many measurement cells that have
% been used for updating the filter.
l = 0;
% Number of partitions of Zk
P = numel(Zk); 
% Allocate some memory
% wp = zeros(1,P);
log_wp = zeros(1,P);
W = zeros(1,P);

% Iterate over the partitions of the measurement set Zk
for p = 1:P
    % Number of cells in partition
    W(p) = numel(Zk(p).P);
    % Allocate memory for dw
    dw = zeros(1,W(p));
    % Iterate over the cells w of the partition p
    for w = 1:W(p)
        % Current set of measurements
        pwZ_k = Zk(p).P(w).W;
        % New measurement cell
        l = l+1;
        % Size of current cell
        absW = size(Zk(p).P(w).W,2);
        
        
        % Iterate over the prediction components
        for j = 1:J_kk1
            % Depending on size of cell, update for extent is performed in
            % different ways.
            if absW < 2
                % Poisson rate
                a_k(1,J_kk1*l+j) = a_kk1(1,j)+absW;
                b_k(1,J_kk1*l+j) = b_kk1(1,j)+1;
                
                % Scalar innovation factor
                S = H_k*P_kk1(:,:,j)*H_k'+1/absW;
                % Gain matrix
                WW = P_kk1(:,:,j)*H_k'*(1/S);
                % Measurement innovation
                innov = pwZ_k-kron(H_k,Id)*m_kk1(:,j);
                % Innovation matrix
                N = (1/S)*(innov*innov');
                % Gaussian mean
                m_k(:,J_kk1*l+j) = m_kk1(:,j)+kron(WW,Id)*innov;
                % Gaussian covariance
                P_k(:,:,J_kk1*l+j) = P_kk1(:,:,j) - WW*S*WW';
                P_k(:,:,J_kk1*l+j) = 0.5*(P_k(:,:,J_kk1*l+j)+P_k(:,:,J_kk1*l+j)');
                % Wishart inverse scale matrix
                V_k(:,:,J_kk1*l+j) = V_kk1(:,:,j)+N;
                V_k(:,:,J_kk1*l+j) = 0.5*(V_k(:,:,J_kk1*l+j)+V_k(:,:,J_kk1*l+j)');
                % Wishart degrees of freedom
                nu_k(1,J_kk1*l+j) = nu_kk1(1,j) + absW;
                
                w_k(1,J_kk1*l+j) = updateGIWweightLog(gam(j),pD(j),...
                    beta_FA,d,absW,S,V_kk1(:,:,j),V_k(:,:,J_kk1*l+j),...
                    nu_kk1(1,j),nu_k(1,J_kk1*l+j))*w_kk1(1,j);
                % Update the weight with Poisson rate likelihood
                w_k(1,J_kk1*l+j) = w_k(1,J_kk1*l+j)*nbinpdf(absW,a_kk1(1,j),b_kk1(1,j)/(b_kk1(1,j)+1));
            else
                % Poisson rate
                a_k(1,J_kk1*l+j) = a_kk1(1,j)+absW;
                b_k(1,J_kk1*l+j) = b_kk1(1,j)+1;
                
                % Centroid measurement
                zbar = mean(pwZ_k,2);
                % Scattering matrix
                zz = pwZ_k-repmat(zbar,1,absW);
                Z = zz*zz';
                
                % Scalar innovation factor
                S = H_k*P_kk1(:,:,j)*H_k'+1/absW;
                % Gain matrix
                WW = P_kk1(:,:,j)*H_k'*(1/S);
                % Centroid measurement innovation
                innov = zbar-kron(H_k,Id)*m_kk1(:,j);
                % Innovation matrix
                N = (1/S)*(innov*innov');
                % Gaussian mean
                m_k(:,J_kk1*l+j) = m_kk1(:,j)+kron(WW,Id)*innov;
                % Gaussian covariance
                P_k(:,:,J_kk1*l+j) = P_kk1(:,:,j) - WW*S*WW';
                P_k(:,:,J_kk1*l+j) = 0.5*(P_k(:,:,J_kk1*l+j)+P_k(:,:,J_kk1*l+j)');
                % Wishart inverse scale matrix
                V_k(:,:,J_kk1*l+j) = V_kk1(:,:,j)+N+Z;
                V_k(:,:,J_kk1*l+j) = 0.5*(V_k(:,:,J_kk1*l+j)+V_k(:,:,J_kk1*l+j)');
                % Wishart degrees of freedom
                nu_k(1,J_kk1*l+j) = nu_kk1(1,j) + absW;
                
                w_k(1,J_kk1*l+j) = updateGIWweightLog(gam(j),pD(j),...
                    beta_FA,d,absW,S,V_kk1(:,:,j),V_k(:,:,J_kk1*l+j),...
                    nu_kk1(1,j),nu_k(1,J_kk1*l+j))*w_kk1(1,j);
                % Update the weight with Poisson rate likelihood
                w_k(1,J_kk1*l+j) = w_k(1,J_kk1*l+j)*nbinpdf(absW,a_kk1(1,j),b_kk1(1,j)/(b_kk1(1,j)+1));
%                 w_k(1,J_kk1*l+j) = w_k(1,J_kk1*l+j)*updateGammaweightLog(a_k(1,J_kk1*l+j),a_kk1(1,j),b_k(1,J_kk1*l+j),b_kk1(1,j));
            end
        end
        % Compute dw
        dw(w) = deltaFunction(absW,1)+sum(w_k(1,J_kk1*l+(1:J_kk1)));
        % Divide the weights by dw
        w_k(1,J_kk1*l+(1:J_kk1)) = w_k(1,J_kk1*l+(1:J_kk1))/dw(w);
    end
    % Replace inf with realmax in dw
    idx = dw==Inf;
    dw(idx) = realmax;
    % Take sum of log(dw) instead of prod(dw) to avoid numerical problems.
    log_wp(p) = sum(log(dw));
end
J_k = J_kk1*(l+1);

% Different normalisation since log(wp) are kept instead of wp.
% See http://en.wikipedia.org/wiki/List_of_logarithmic_identities
if P == 0
    wp = [];
elseif P == 1
    wp = 1;
else
    log_wp = log_wp-(log_wp(1)+log(1+sum(exp(log_wp(2:end)-log_wp(1)))));
    wp = exp(log_wp);
end

% figure(111),clf
% plot(wp,'.')

prev = J_kk1;
for p = 1:length(wp)
    w_k(1,prev+(1:J_kk1*W(p))) = w_k(1,prev+(1:J_kk1*W(p)))*wp(p);  
    prev = prev + J_kk1*W(p);
end

w_k = w_k(1,1:J_k);
a_k = a_k(1,1:J_k);
b_k = b_k(1,1:J_k);
m_k = m_k(:,1:J_k);
P_k = P_k(:,:,1:J_k);
nu_k = nu_k(1,1:J_k);
V_k = V_k(:,:,1:J_k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [df] = deltaFunction(x,y)
% this is just an implementation of Diracs delta...
if x == y
    df = 1;
else
    df = 0;
end 
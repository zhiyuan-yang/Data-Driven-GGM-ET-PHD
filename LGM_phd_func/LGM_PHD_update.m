function [w_k,a_k,b_k,m_k,P_k,J_k] = ...
    LGM_PHD_update(w_kk1,a_kk1,b_kk1,m_kk1,P_kk1,J_kk1,Zk,model,filter,ego,GMModel)
% update process of PHD Filter
nx = size (m_kk1,1); %[x y vx vy l w] 6

% UT parameter
alpha = filter.LGM_phd.alpha;
beta = filter.LGM_phd.beta;
kappa = filter.LGM_phd.kappa;
% measurement noise covariance
R = model.obs.R; 
% Probabilities and detection
p_D_k = model.obs.pD;

% Clutter rate
beta_FA = model.obs.lambda_c * model.obs.pdf_c;

% Allocate some memory
nZk = size(Zk,2);
w_k = zeros(1,(nZk+1)*J_kk1);
a_k = zeros(1,(nZk+1)*J_kk1);
b_k = zeros(1,(nZk+1)*J_kk1);
m_k = zeros(nx,(nZk+1)*J_kk1);
P_k = zeros(nx,nx,(nZk+1)*J_kk1);

% The expected number of generated measurements
gam(1:J_kk1) = a_kk1(1:J_kk1)./b_kk1(1:J_kk1);   %beta_D;

% Allocate memory for probability of detection
pD = p_D_k*ones(1,J_kk1);
for j = 1:J_kk1
    % Gaussian component corresponding to no detections
    w_k(1,j) = (1-(1-exp(-gam(j)))*pD(j))*w_kk1(1,j);
    a_k(1,j) = a_kk1(1,j);
    b_k(1,j) = b_kk1(1,j);
    m_k(:,j) = m_kk1(:,j);
    P_k(:,:,j) = P_kk1(:,:,j);
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
            % likelihood function
            PL = zeros(1,absW);
            
            a_k(1,J_kk1*l+j) = a_kk1(1,j)+absW;
            b_k(1,J_kk1*l+j) = b_kk1(1,j)+1;
            mx = m_kk1(:,j);
            Px= P_kk1(:,:,j);
            for nz =1:absW          %Iterate over the measurements one by one
        
                [xPts,w_m,w_c] = UT(mx,Px,alpha,beta,kappa); 
                ego_pos = ego.position(1:2);
                %note that phi should be calculated in sensor-coordinate
                phi = atan2(xPts(2,:)-ego_pos(2,:),xPts(1,:)-ego_pos(1,:));
                zPts = f_meas_pseudo_LGM(xPts, pwZ_k(:,nz), phi, GMModel); %zPts 2*13 vector
                nPts = size(zPts, 2);
                
                zhat_nz = [0;0];
                for i = 1:nPts
                    zhat_nz = zhat_nz + (w_m(i) * zPts(:,i));   %zhat_nz 2*1 vector
                end
                
                S_nz = zeros(2,2);
                P_xz = zeros(nx,1);
                for i = 1:nPts
                    S_nz = S_nz + w_c(i) * ((zPts(:,i) - zhat_nz) * (zPts(:,i) - zhat_nz)'); %2*2 matrix
                    P_xz = P_xz + w_c(i) * ((xPts(:,i) - mx) * (zPts(:,i) - zhat_nz)');     %6*2 matrix
                end
                S_nz = S_nz + R;
                % kalman gain
                K = P_xz / S_nz;
                
                mx = mx + K * (pwZ_k(:,nz) - zhat_nz);
                Px = Px - K * S_nz * K';
                PL(nz) = mvnpdf(pwZ_k(:,nz)',zhat_nz',S_nz);
            end
            m_k(:,J_kk1*l+j) = mx;
            P_k(:,:,J_kk1*l+j) = Px;
            P_k(:,:,J_kk1*l+j) = 0.5*(P_k(:,:,J_kk1*l+j)+P_k(:,:,J_kk1*l+j)');
            w_k(1,J_kk1*l+j) = exp(-gam(j)+log(pD(j))+absW*log(gam(j))+sum(log(PL))-absW*log(beta_FA))*w_kk1(1,j);   
            w_k(1,J_kk1*l+j) = w_k(1,J_kk1*l+j)*updateGammaweightLog(a_k(1,J_kk1*l+j),a_kk1(1,j),b_k(1,J_kk1*l+j),b_kk1(1,j));
%             w_k(1,J_kk1*l+j) = w_k(1,J_kk1*l+j)*nbinpdf(absW,a_kk1(1,j),b_kk1(1,j)/(b_kk1(1,j)+1));
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
test = length(wp);
for p = 1:1:test
    w_k(1,prev+(1:J_kk1*W(p))) = w_k(1,prev+(1:J_kk1*W(p)))*wp(p);  
    prev = prev + J_kk1*W(p);
end

w_k = w_k(1,1:J_k);
a_k = a_k(1,1:J_k);
b_k = b_k(1,1:J_k);
m_k = m_k(:,1:J_k);
P_k = P_k(:,:,1:J_k);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [df] = deltaFunction(x,y)
% this is just an implementation of Diracs delta...
if x == y
    df = 1;
else
    df = 0;
end



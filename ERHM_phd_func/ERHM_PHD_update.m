function [w_k,a_k,b_k,m_k,P_k,J_k] = ...
    ERHM_PHD_update(w_kk1,a_kk1,b_kk1,m_kk1,P_kk1,J_kk1,Zk,model,filter)

x_dim = model.RHM.x_dim; % kin.x_dim =4
z_dim = model.obs.z_dim;
nx = size (m_kk1,1); %kin.x_dim(4)+l_dim(3)=7
% UT parameter
alpha = filter.ERHM_phd.alpha;
beta = filter.ERHM_phd.beta;
kappa = filter.ERHM_phd.kappa;
% measurement noise covariance
R = model.obs.R; 
% Probabilities and detection
p_D_k = model.obs.pD;

% Clutter rate
beta_FA = model.obs.lambda_c * model.obs.pdf_c;

% scale factor
s_mean = model.RHM.s_mean;
s_var = model.RHM.s_var;

% Allocate some memory
nZk = size(Zk,2);
w_k = zeros(1,(nZk+1)*J_kk1);
a_k = zeros(1,(nZk+1)*J_kk1);
b_k = zeros(1,(nZk+1)*J_kk1);
m_k = zeros(size(m_kk1,1),(nZk+1)*J_kk1);
P_k = zeros(size(P_kk1,1),size(P_kk1,2),(nZk+1)*J_kk1);

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
            m_center = m_kk1(1:2,j);
            for nz =1:absW          %Iterate over the measurements one by one
                %增补后的状态
                mx_aug = [mx',s_mean,zeros(1,z_dim)]';
                nx_aug = length(mx_aug); %kin.x_dim(4)+l_dim(3)+s(1)+noise(2)=10
                Px_aug = blkdiag(Px,s_var,R);
                theta = atan2(pwZ_k(2,nz) - m_center(2), pwZ_k(1,nz) - m_center(1));
                
                [xPts,w_m,w_c] = UT(mx_aug,Px_aug,alpha,beta,kappa);
                zPts = f_meas_pseudo_ellipse(xPts(1:nx,:), xPts(nx+1:nx_aug,:), pwZ_k(:,nz),theta, x_dim);
                nPts = size(zPts, 2);
                
                zhat_nz = 0;
                for i = 1:nPts
                    zhat_nz = zhat_nz + (w_m(i) * zPts(i));
                end
                
                S_nz = 0;
                P_xz = zeros(nx_aug,1);
                for i = 1:nPts
                    S_nz = S_nz + w_c(i) * ((zPts(i) - zhat_nz) * (zPts(i) - zhat_nz)');
                    P_xz = P_xz + w_c(i) * ((xPts(:,i) - mx_aug) * (zPts(i) - zhat_nz)');
                end
                % kalman gain
                K = P_xz / S_nz;
                
                mx_aug = mx_aug + K * (0 - zhat_nz);
                Px_aug = Px_aug - K * S_nz * K';
                PL(nz) = normpdf(0,zhat_nz,S_nz);
                
                mx = mx_aug(1:nx);
                Px = Px_aug(1:nx,1:nx);
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
for p = 1:length(wp)
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

function pseudoMeasurement = f_meas_pseudo_ellipse(x, noise, z , theta, x_dim)
numberOfSigmaPoints = size(x, 2);
pseudoMeasurement = zeros(1, numberOfSigmaPoints);
for j = 1: numberOfSigmaPoints
    s = noise(1,j);                                 %scale factor
    v = noise(2:3,j);                               %measurement noise
    m = x(1:2,j);                                   %the center position
    L = [x(x_dim+1,j),0;
        x(x_dim+3,j),x(x_dim+2,j)];
    Ellipara = calculateEllipara(L*L');
    
%     theta = atan2(z(2) - m(2), z(1) - m(1))+2*pi;
    R = calcPolarFunc(theta, Ellipara);
    e = [cos(theta); sin(theta)];
    pseudoMeasurement(j) = (s * R)^2 + 2 * s * R * e' * v + norm(v)^2 - (norm( m - z ))^2 ;
end

function R = calcPolarFunc(theta,Ellipara)
a = Ellipara(1);
b = Ellipara(2);
phi = Ellipara(3);

R = a*b/sqrt(a*sin(theta-phi).^2+b*cos(theta-phi).^2);
R = real(R);


function Ellipara = calculateEllipara(A)
E = eig(A);
%椭圆的长短轴a,b分别沿着矩阵A的两个特征向量的方向，
%而两个与之对应的特征值分别是半长轴和半短轴的长度的平方的倒数
a = max(sqrt(E));
b = min(sqrt(E));
% if a-b < 1e-5
%     phi =0;
% else
%     phi = real(asin(2*A(1,2)/(a^2-b^2)))/2+2*pi;
% end
phi = 0;
Ellipara(1)=a;
Ellipara(2)=b;
Ellipara(3)=phi;

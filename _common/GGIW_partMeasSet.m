function [Zp,dG] = GGIW_partMeasSet(Z,w,m,a,b,nu,V,model,filter)

max_distance = model.max_distance;
min_distance = model.min_distance;
F = model.GGIW.F_k;

plotPartitions = filter.plot_partitions;

if isempty(Z)
    Zp = [];
elseif size(Z,2) == 1
    Zp(1).P(1).W = Z;
    dG = 0;
else
    [Zp,dG] = partition(Z,max_distance,min_distance);
    
    P_W = [Zp.N_W];
%     [Zp,dG,P_W] = subPartitionKmeans(Zp,dG,a,b,m,w);

    Np = numel(Zp);
    P_pred = predictionPartition(model,Z,w,m,nu,V,F);
    if ~isempty(P_pred)
        Zp(Np+1).P = P_pred;
        Zp(Np+1).N_W = numel(P_pred);
        P_W(Np+1) = numel(P_pred);
        [Npp] = findEqualPartitions(Zp,P_W,Np+1);
        if Npp == Np+1
            Np = Npp;
            dG(end+1) = -1;
        else
            Zp = Zp(1:Np);
        end
    end
%     P_EM = EMpartition(model,Z,w,m,a,b,nu,V,F);
%     if ~isempty(P_EM)
%         Zp(Np+1).P = P_EM;
%         Zp(Np+1).N_W = numel(P_EM);
%         P_W(Np+1) = numel(P_EM);
%         [Npp] = findEqualPartitions(Zp,P_W,Np+1);
%         if Npp == Np+1
%             
%             dG(end+1) = -2;
%         else
%             Zp = Zp(1:Np);
%         end
%     end
end

if ~isunix && plotPartitions && ~isempty(Z)
%     global xMin xMax yMin yMax
    Ncols = 4;
    figure(67),clf,%whitebg('k')
    subplot(ceil((numel(Zp)+1)/Ncols),Ncols,1)
    plot(Z(1,:),Z(2,:),'o','markerfacecolor','b','markeredgecolor','k')
    axis equal
    title(['Measurements ' num2str(length(Z(1,:)))])
    xlabel('X [m]')
    ylabel('Y [m]')
%     axis([xMin xMax yMin yMax])
    for p = 1:numel(Zp)
        W = numel(Zp(p).P);
        col = jet(W);
        subplot(ceil((numel(Zp)+1)/Ncols),Ncols,p+1)
        for w = 1:W
            plot(Zp(p).P(w).W(1,:),Zp(p).P(w).W(2,:),'o','markerfacecolor',col(w,:),'markeredgecolor','k')
            hold on
            axis equal
            title(['d = ' num2str(dG(p))])
            xlabel('X [m]')
            ylabel('Y [m]')
            Xmax = max(Zp(p).P(w).W,[],2);
            Xmin = min(Zp(p).P(w).W,[],2);
            e = 1;
            plot([Xmin(1)-e Xmin(1)-e Xmax(1)+e Xmax(1)+e Xmin(1)-e],...
                [Xmin(2)-e Xmax(2)+e Xmax(2)+e Xmin(2)-e Xmin(2)-e],'k')
            text(0.5*(Xmin(1)+Xmax(1)),0.5*(Xmin(2)+Xmax(2)),num2str(length(Zp(p).P(w).W(1,:))))
%             axis([xMin xMax yMin yMax])
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Zp,dG] = partition(Z,max_distance,min_distance)
% Takes a scan L and divides it into subgroups

% Number of points in scan
N = size(Z,2);
% X and Y positions
X = Z(1,:); Y = Z(2,:);
% Matrix with point to point distances
distMat = sqrt((repmat(X,N,1)-repmat(X',1,N)).^2+...
    (repmat(Y,N,1)-repmat(Y',1,N)).^2);
% Get sorted list with point to point distances
[dMsort_tmp, idxSort_tmp] = sort(distMat(:));
% Get unique distances
[dMsort,uniIdx,tmp] = unique(dMsort_tmp);
% Get unique indexes
idxSort = idxSort_tmp(uniIdx);
% Find idx to distance larger than max_distance
dMidx = find(dMsort > max_distance,1,'first');
% Extract the distances
if isempty(dMidx)
    dG = dMsort(1:end);
    dMidx = length(dG);
else
    dG = dMsort(1:dMidx);
end
% Obtain the point to point correnspondences of idxSort
[I,J] = ind2sub(size(distMat),idxSort);
% Extract the short ones
I = I(1:dMidx); J = J(1:dMidx);
% Allocate memory for which points are in the same cluster
cluster = zeros(size(distMat));
% Allocate memory
Cmat = zeros([size(distMat,1) dMidx-1]);
cellLabel = 1:N;
Zp = struct('P',repmat({[]},1,N));
partition = 1;
cells = struct('W',repmat({[]},1,N));
for iz=1:N
    cells(iz).W = Z(:,iz);
    cells(iz).idx = iz;
end
% Set the first cluster corresponding to one meas in each cell
Cmat(:,1) = 1:size(distMat,1);
Zp(partition).P = cells;
Zp(partition).N_W = N;

if dMidx == 2
    if dG(2) < min_distance
        dG = min_distance;
        Zp(1).P = struct('W',Z);
        Zp(1).N_W = 1;
        Zp(1).P(1).idx = 1:size(Z,2);
        Zp = Zp(1);
    else
        dG = min(dG(2),max_distance);
        Zp = Zp(1);
    end
else
    distG = 0;
    dG = (dG(1:end-1)+dG(2:end))/2;
    for k = 2:dMidx-1
        % Make sure both points not already included
        if ~cluster(I(k),J(k)) % if one, then already in the same
            % Increment partition counter
            partition = partition+1;
            % Set the new cluster to the previous one
            Cmat(:,partition) = Cmat(:,partition-1);
            % Save this distance
            distG = [distG dG(k)]; %(dG(k)+dG(k+1))/2];
            % Index to points that I and J have previously been clustered to
            Iprev = find(cluster(I(k),:));
            Jprev = find(cluster(J(k),:));
            % Set those points to one
            idx1 = unique([I(k) J(k) Iprev Jprev]);
            cluster(idx1,idx1) = 1;
            % Fing group labels
            grpLabel = sort(Cmat([I(k) J(k)],partition));
            % Corresponding indexes
            gi1 = find(cellLabel==grpLabel(1));
            gi2 = find(cellLabel==grpLabel(2));
            % Set the cell labels
            cellLabel = setdiff(cellLabel,grpLabel(2));
            % Find minimal group label
            Cmidx = min(grpLabel);
            % Set group label
            Cmat(idx1,partition) = Cmidx;
            
            % Put the corresponding measurements in the same cell
            cells(gi1).W = [cells(gi1).W cells(gi2).W];
            cells(gi1).idx = [cells(gi1).idx cells(gi2).idx];
            % Remove the second cell
            cells = cells(setdiff(1:numel(cells),gi2));
            
            Zp(partition).P = cells;
            Zp(partition).N_W = numel(cells);
        end
    end
    dG = distG;
    % Take the distances larger than min_distance
    idx2 = find(dG>min_distance);
    if isempty(idx2)
        dG = min_distance;
        Zp = Zp(partition);
    else %if idx2(1)>1
        dG = dG(idx2);
        dG = [min_distance dG];
        idx2 = [idx2(1)-1 idx2];
        % set the partitions
        Zp = Zp(idx2);
    end
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Zp,dG,P_W] = subPartitionKmeans(Zp,dG,a,b,m,w)
% % Poisson rate for #meas per target, and maximum number of targets
% global beta_D Nx_max d
% % Number of partitions
% N_P = numel(Zp);
% % Counter over number of new partitions
% ip_new = N_P;
% % Number of cells per partition
% P_W = [Zp.N_W];
% 
% idx_w = find(w>0.5);
% if ~isempty(idx_w)
%     a=a(idx_w);
%     b=b(idx_w);
%     m=m(:,idx_w);
% end
% 
% % Iterate over the partitions
% for ip = 1:N_P
%     % Current partition
%     Zp_P = Zp(ip).P;
%     % Number of cells in the partition
%     N_W = numel(Zp_P);
%     P_W(ip) = N_W;
%     % Iterate over the cells
%     for iw = 1:N_W
%         % Current cell
%         Zp_P_W = Zp_P(iw).W;
%         Zp_P_idx = Zp_P(iw).idx;
%         % Size of current cell
%         absW = size(Zp_P_W,2);
%         % Cell center point
%         Wc = mean(Zp_P_W,2);
%         % Find closest component
%         if ~isempty(idx_w)
%             [mydummy,idxc]=min((m(1,:)-Wc(1)).^2+(m(2,:)-Wc(2)).^2);
%             gamma_hat = max(3,a(idxc)/b(idxc));
%             Nx_hat = round(absW/(gamma_hat));
%         else
%             gamma_hat = mean(max(3,a./b));
%             Nx_hat = round(absW/gamma_hat);
%         end
%         
%         % Nx_hat can not be larger than absW
%         Nx_hat = min(Nx_hat,absW);
%         % Create new partition if more than one
%         if Nx_hat > 1
%             % Increase partition counter
%             ip_new = ip_new+1;
%             % Copy the partition
%             Zp(ip_new).P = Zp_P;
%             % Split the measurements using K-means
%             idx = kmeans(Zp_P_W',Nx_hat,'Replicates',10);
%             % Index to cells that are to be changed or added
%             Widx = [iw N_W-1+(2:Nx_hat)];
%             % Reset the cell that is being split and add new cells
%             for iw2 = 1:length(Widx)
%                 Zp(ip_new).P(Widx(iw2)).W = Zp_P_W(:,idx==iw2);
%                 Zp(ip_new).P(Widx(iw2)).idx = Zp_P_idx(:,idx==iw2);
%             end
%             P_W(ip_new) = numel(Zp(ip_new).P);
%             Zp(ip_new).N_W = P_W(ip_new);
%             % add the distance
%             dG = [dG dG(ip)];
%         end
%         [ip_new] = findEqualPartitions(Zp,P_W,ip_new);
%     end
% end
% Zp = Zp(1:ip_new);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ip_new] = findEqualPartitions(Zp,P_W,ip_new)
same = 0;
for ip2 = find(P_W(1:ip_new-1)==P_W(ip_new))
    for iw_new = 1:P_W(ip_new)
        for iw_old=1:P_W(ip2)
            same = same+isempty(setxor(Zp(ip2).P(iw_old).idx,Zp(ip_new).P(iw_new).idx));
        end
    end
end
% If the number of same cells are equal to number of cells,
% then the partitions are equal.
if same==P_W(ip_new)
    ip_new=ip_new-1;
%     display('same')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P_pred = predictionPartition(model,Z,w,m,nu,V,F)

% Function that creates a partition based on the predictions of extracted
% components.

% Sampling time and temporal decay
Ts = model.GGIW.Ts;
tau = model.GGIW.tau;

% Smallest number of Wishart degrees of freedom allowed
nu_min = model.GGIW.v_min;

% Number of measurements
Nz = size(Z,2);
Iz = 1:Nz;
% Allocate memory
P_pred = struct('W',[],'idx',[]);
% Size of extension
d = size(V,1);
Id = eye(d);
% Sort weights and find the extracted ones
[w,iw] = sort(w,'descend');
IK = iw(w>0.5);
% Counter for the cells
iw = 0;

% Prediction Partition gate
predPartGate=chi2inv(0.99,d);

for ik = IK
    % Predict and find extension estimate
    mp = kron(F,Id)*m(:,ik);
    nup = max(exp(-Ts/tau)*nu(ik),nu_min);
    Vp = ((nup-d-1)/(nu(ik)-d-1))*V(:,:,ik);
    Ext = Vp/max(1,nup-2*d-2);
    % Find measurements inside extension
    Zm = Z(:,Iz)-repmat(mp(1:2)',1,length(Iz));
    idx = find(diag(Zm'*(Ext\Zm))' < predPartGate);
    % Add new cell
    if ~isempty(idx)
        iw = iw+1;
        P_pred(iw).W = Z(:,Iz(idx));
        P_pred(iw).idx = Iz(idx);
    end
    % Remove cell from set of measurements
    Iz = setdiff(Iz,Iz(idx));
end

if iw == 0
    P_pred = [];
elseif ~isempty(Iz)
    % Put remaining measurements in single cells
    for ik = 1:length(Iz)
        iw = iw+1;
        P_pred(iw).W = Z(:,Iz(ik));
        P_pred(iw).idx = Iz(ik);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P_EM = EMpartition(model,Zk,w,m,a,b,nu,V,F)

% Function that computes a partition by performing EM for Gaussian Mixtures
% (mean and extension) to the measurements.

% Surveilance area
xMin=model.obs.xrange(1);
xMax=model.obs.xrange(2);
yMin=model.obs.yrange(1);
yMax=model.obs.yrange(2);
% Sampling time, temporal decay, prediction for Poisson rates
Ts = model.GGIW.Ts;
tau = model.GGIW.tau;
eta_k =model.GGIW.eta_k;

% Smallest number of Wishart degrees of freedom allowed
nu_min = model.GGIW.v_min;
% Size of extension
d = size(V,1);
Id = eye(d);
% Sort weights and find the extracted ones
[w,iw] = sort(w,'descend');
IK = iw(w>0.5);
% Number of measurements
Nz = size(Zk,2);

N = length(IK);
% Initialise means, covariances, mixing coefficients
mu_k = zeros(d,N);
Sig_k = zeros(d,d,N);
pi_k = zeros(1,d);
for k = 1:N
    % Predict and find extension estimate
    mp = kron(F,Id)*m(:,IK(k));
    mu_k(:,k) = mp(1:2);
    nup = max(exp(-Ts/tau)*nu(IK(k)),nu_min);
    Vp = ((nup-d-1)/(nu(IK(k))-d-1))*V(:,:,IK(k));
    Sig_k(:,:,k) = Vp/max(1,nup-2*d-2);
    
    ap = a(k)/eta_k;
    bp = b(k)/eta_k;
    pi_k(k) = ap/bp;
    % Use the degrees of freedom as estimate instead...
    % compGammaFromExt(Sig_k(:,:,k));
end
% Add additional Gaussian component supposed to catch clutter
% Mean position(s)
perDim = 1;
[Xe,Ye] = meshgrid(xMin+0.5*(xMax-xMin)/perDim:(xMax-xMin)/perDim:xMax-0.5*(xMax-xMin)/perDim,...
    yMin+0.5*(yMax-yMin)/perDim:(yMax-yMin)/perDim:yMax-0.5*(yMax-yMin)/perDim);
% Number of components to add
Ne = perDim^2;
for k = 1:Ne
    mu_k(:,N+k) = [Xe(k);Ye(k)];
    Sig_k(1:2,1:2,N+k) = diag([0.5*(xMax-xMin)/(3*perDim) 0.5*(xMax-xMin)/(3*perDim)].^2);
    pi_k(N+k) = 1e-9;
end
% Total number of components
nc = N+Ne;
% Normalise mixing coefficients
pi_k = pi_k/sum(pi_k);
% Allocate memory for EM-algorithm
GAM = zeros(Nz,nc);
lik = zeros(Nz,nc);
% Maximum number of iterations
niter = 20;
% Convergence criterion for log-likelihood
convCrit = 0.01;
Lold = -inf;
for k = 1:niter
    % E-step
    for ic = 1:nc
        GAM(:,ic) = pi_k(ic)*mvnpdf(Zk',mu_k(:,ic)',Sig_k(:,:,ic));
    end
    GAM(GAM==inf) = 1e100; % To avoid infinite probabilities
    GAM=GAM./repmat(sum(GAM,2),1,nc);

    % M-step
    N_k = sum(GAM,1);
    N_k = max(N_k,1e-9); % To avoid N_k=0
    for ic = 1:nc
        mu_k(:,ic) = sum(Zk.*repmat(GAM(:,ic)',2,1),2)/N_k(ic);
        Z_mu = (Zk-repmat(mu_k(:,ic),1,Nz)).*repmat(sqrt(GAM(:,ic)'),2,1);
        Sig_k(:,:,ic) = Z_mu*Z_mu'/N_k(ic);
        Sig_k(:,:,ic) = 0.5*(Sig_k(:,:,ic)+Sig_k(:,:,ic)');
        
        [eig_v,eig_d]=eig(Sig_k(:,:,ic));
        if min(diag(eig_d)) < 1e-3
            % Matrix is nearly singular, perform Haxxx
            Sig_k(:,:,ic) = eig_v*diag(max(diag(eig_d),1e-3))*eig_v';
        end
    end
    pi_k = N_k/Nz;
    % Evaluate log likelihood
    for ic = 1:nc
        lik(:,ic) = pi_k(ic)*mvnpdf(Zk',mu_k(:,ic)',Sig_k(:,:,ic));
    end
    Lnew= sum(log(sum(lik,2)));
    
    if k>1 && abs(Lnew-Lold)<convCrit
        break
    else
        Lold = Lnew;
    end
end
% Find associations
[mydummy,idx]=max(GAM,[],2);
% if all are associated to clutter
if sum(idx==size(GAM,2))==Nz
    P_EM=[];
else
    % Allocate memory
    P_EM = struct('W',{[]},'idx',{[]});
    iw = 0;
    % Cells corresponding to targets
    for kk = 1:size(GAM,2)-1
        if sum(idx==kk)>0 % Make sure that empty cells are not created
            iw = iw+1;
            P_EM(iw).W = Zk(:,idx==kk);
            P_EM(iw).idx = find(idx==kk);
        end
    end
    cidx = find(idx==size(GAM,2));
    Nc = sum(idx==size(GAM,2));
    for kk = 1:Nc
        iw = iw+1;
        P_EM(iw).W = Zk(:,cidx(kk));
        P_EM(iw).idx = cidx(kk)';
    end
end
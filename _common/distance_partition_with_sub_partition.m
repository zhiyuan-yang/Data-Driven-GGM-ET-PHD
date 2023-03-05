function [Zp,dG] =distance_partition_with_sub_partition(Z,model,filter)

max_distance = model.max_distance;
min_distance = model.min_distance;
lambda = model.lambda;
N_low = filter.N_low_FA;

if isempty(Z)
    Zp = [];
elseif size(Z,2) == 1
    Zp(1).P(1).W = Z;
    Zp(1).N_W=1;
    dG = 0;
else
    [Zp,dG] = partition(Z,max_distance,min_distance);
    [Zp,dG] = split(Zp,dG,lambda);
    [Zp,dG] = mergePartitions(Z,Zp,dG,N_low);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Zp,dG] = partition(Z,max_distance,min_distance)
% Partitions a set of measurements

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
clust = zeros(size(distMat));
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
        if ~clust(I(k),J(k)) % if one, then already in the same
            % Increment partition counter
            partition = partition+1;
            % Set the new cluster to the previous one
            Cmat(:,partition) = Cmat(:,partition-1);
            % Save this distance
            distG = [distG dG(k)]; %(dG(k)+dG(k+1))/2];
            % Index to points that I and J have previously been clustered to
            Iprev = find(clust(I(k),:));
            Jprev = find(clust(J(k),:));
            % Set those points to one
            idx1 = unique([I(k) J(k) Iprev Jprev]);
            clust(idx1,idx1) = 1;
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


function [Zp,dG] = mergePartitions(Z,Zp,dG,N_low)

% Number of partitions
Np = numel(Zp);
% counter for new partitions
ip_new = Np;
% Iterate over all partitions
for p = 1:Np
    % Current partition
    Z_P = Zp(p).P;
    % Number of cells in current partition
    Nw = numel(Z_P);
    
    for nlow = 1:N_low
        ip_new = ip_new+1;
        % Allocate memory
        idx_FA = [];
        idx_T = [];
        % Iterate over cells
        for w = 1:Nw
            if size(Zp(p).P(w).W,2)<=nlow
                % If there is less than N_low measurements in the cell
                idx_FA = [idx_FA Zp(p).P(w).idx];
            else
                % If there are more than N_low measurement
                idx_T = [idx_T w];
            end
        end
        % Keep the cells with more than N_low measurements
        Z_P_new = Z_P(idx_T);
        if ~isempty(idx_FA)
            % Merge all the other cells
            Z_P_new(length(idx_T)+1).W = Z(:,idx_FA);
            Z_P_new(length(idx_T)+1).idx = idx_FA;
        end
        % Save the new partition
        Zp(ip_new).P = Z_P_new;
        Zp(ip_new).N_W = numel(Z_P_new);
        dG(ip_new) = dG(p);
        % Check for equal partitions
        [ip_new] = findEqualPartitions(Zp,[Zp.N_W],ip_new);
    end
end
Zp = Zp(1:ip_new);
dG = dG(1:ip_new);

function [Zp,dG] = split(Zp, dG, lambda)
num_p = length(Zp);
ip_new = num_p + 1;
dele_idx = [];  %index of which partition should be deleted
for i = 1:1:num_p
    curr_p = Zp(i);
    num_w = curr_p.N_W;
    flag_split = 0;
    for j = 1:1:num_w
        idx_newp = 1;
        num_meas = length(curr_p.P(j).idx);
        [~ , max_n] = max(poisspdf(num_meas,lambda*[1:5]));
        if max_n > 1  %cell need to be splitted  
            temp_meas = curr_p.P(j).W;
            temp_idx = curr_p.P(j).idx;
            idx = kmeans(temp_meas',max_n);
            if flag_split == 0  %no previous cell is splitted
                if j>1  %there is a previou cell
                    for k = 1:j-1
                        Zp(ip_new).P(idx_newp) = curr_p.P(k); %copy all the previous measurement cell
                        idx_newp = idx_newp + 1;
                    end
                end
                for k = 1:max_n
                    Zp(ip_new).P(idx_newp).W = temp_meas(:,find(idx == k));
                    Zp(ip_new).P(idx_newp).idx = temp_idx(find(idx == k));  %split measurement cell
                    idx_newp = idx_newp + 1;
                end
                if j<num_w  %there is a followup cell
                    for k = j+1:num_w
                        Zp(ip_new).P(idx_newp) = curr_p.P(k);%copy all the measurement cell below
                        idx_newp = idx_newp + 1;
                    end
                end
                Zp(ip_new).N_W = numel(Zp(ip_new).P);
                dG(ip_new) = dG(i);
            else  %previous cell has been splitted
                num_w_newp = length(Zp(ip_new).P);
                idx_newp = num_w_newp - num_w + j;
                for k = 1:max_n 
                    Zp(ip_new).P(idx_newp).W = temp_meas(:,find(idx == k));
                    Zp(ip_new).P(idx_newp).idx = temp_idx(find(idx == k));
                    idx_newp = idx_newp + 1;
                end
                if j<num_w 
                    for k = j+1:num_w
                        Zp(ip_new).P(idx_newp) = curr_p.P(k);
                        idx_newp = idx_newp + 1;
                    end
                end
                Zp(ip_new).N_W = numel(Zp(ip_new).P);
                dG(ip_new) = dG(i);
            end
            flag_split = 1;
        end
    end
    if flag_split == 1
        [temp] = findEqualPartitions(Zp,[Zp.N_W],ip_new);  %find if new partition is redundent
        if temp ~= ip_new
            Zp(ip_new) = [];  %new partition is redundent delete it
        else
            ip_new = ip_new + 1;
        end
        dele_idx = [dele_idx,i];  %partition with split measurement 
    end
end
if ~isempty(dele_idx)
    Zp(dele_idx) = [];
end


function [ip_new] = findEqualPartitions(Zp,P_W,ip_new)
% Find partitions with same number of cells
idx_same_Nw = find(P_W(1:ip_new-1)==P_W(ip_new));
Nsame = length(idx_same_Nw);
same = zeros(1,Nsame);
for i_same = 1:Nsame
    ip2 = idx_same_Nw(i_same);
    for iw_new = 1:P_W(ip_new)
        for iw_old=1:P_W(ip2)
            same(i_same) = same(i_same)+isempty(setxor(Zp(ip2).P(iw_old).idx,Zp(ip_new).P(iw_new).idx));
        end
    end
end
% If the number of same cells are equal to number of cells,
% then the partitions are equal.
if sum(same==P_W(ip_new))>0
    ip_new=ip_new-1;
%     display('same')
end


function [wb,ab,bb,mb,Pb,Jb] = LGM_reduction(w,a,b,m,P,J,filter)

% Function that prunes GGM components with weights w<T and merges
% components whose pairwise KL-difference is <U. At most Jmax components
% are returned.
T = filter.LGM_phd.T;
U = filter.LGM_phd.U;
Jmax = filter.LGM_phd.J_max;

state_dim = size(m,1);      

% Allocate memory
wb=zeros(1,J);
ab=zeros(1,J);
bb=zeros(1,J);
mb=zeros(state_dim,J);
Pb=zeros(state_dim,state_dim,J);


% Do not use components with too low weights
I = find(w > T);
el = 0;
while ~isempty(I)    % 较大的权重
    el = el+1;
    % Find maximum weight
    [~,jtmp] = max(w(I));
    j = I(jtmp);      % 找到权重最大的项
    
    iPj = inv(P(1:2,1:2,j));
    L = [];
    % Iterate over components
    for i = I
        val = (m(1:2,i)-m(1:2,j))'*iPj*(m(1:2,i)-m(1:2,j));
        if val <= U
            L = [L i];
        end
    end
    
    % Save the merged weight
    wb(el) = sum(w(L));
    % Merge the components
    ab(el) = sum(w(L).*a(L))/wb(el);
    bb(el) = sum(w(L).*b(L))/wb(el);
    mb(:,el) = wsumvec(w(L),m(:,L),state_dim)/wb(el);
    Pb(:,:,el) = wsummat(w(L),P(:,:,L),state_dim)/wb(el);  
    
    % Remove the merged components from the set of components
    I = setdiff(I,L);
end

% To get the same sum of weights
sumwb = sum(wb);
% Check if there are too many elements
Jb = min(el,Jmax);
% Sort the merged weights
[~,i] = sort(wb,'descend');
i = i(1:Jb);
% Keep components ordered after weight
wb = wb(1,i);
wb = sumwb*wb/sum(wb);
ab = ab(1,i);
bb = bb(1,i);
mb = mb(:,i);
Pb = Pb(:,:,i);

function out = wsumvec(w,vecstack,xdim)
    wmat = repmat(w,[xdim,1]);
    out  = sum(wmat.*vecstack,2);

function out = wsummat(w,matstack,xdim)
    w = reshape(w,[1,1,length(w)]);
    wmat = repmat(w,[xdim,xdim,1]);
    out = sum(wmat.*matstack,3);
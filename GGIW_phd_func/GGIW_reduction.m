function [wb,ab,bb,mb,Pb,vb,Vb,Jb] = GGIW_reduction(w,a,b,m,P,v,V,J,model,filter)

% Function that prunes GIW components with weights w<T and merges
% components whose pairwise KL-difference is <U. At most Jmax components
% are returned.
T = filter.GGIW_phd.T;
U_G = filter.GGIW_phd.U_G;
U_N = filter.GGIW_phd.U_N;
U_IW = filter.GGIW_phd.U_IW;
Jmax = filter.GGIW_phd.J_max;
plotReduction = filter.GGIW_phd.plot_reduction;

d   = size(V,1);      % 2*2*N
n_x = size(m,1);      % 目标状态维数
s   = n_x/d;

% Allocate memory
wb=zeros(1,J);
ab=zeros(1,J);
bb=zeros(1,J);
mb=zeros(n_x,J);
Pb=zeros(s,s,J);
vb=zeros(1,J);
Vb=zeros(d,d,J);

% Do not use components with too low weights
I = find(w>T);
fig_ctr = 3;
l = 0;
while ~isempty(I)    % 较大的权重
    l = l+1;
    % Find maximum weight
    [~,jtmp] = max(w(I));
    j = I(jtmp);      % 找到权重最大的项
    
    if plotReduction
        fig_ctr=fig_ctr+1;
        figure(fig_ctr),clf,hold on,axis equal,%axis([-20 20 -5 20])
        plotGGIWcomponent(a(j),b(j),m(:,j),P(:,:,j),v(j),V(:,:,j),'b',6)
    end
    
    Pj = kron(P(:,:,j),V(:,:,j))/(v(j)+s-s*d-2);
    L = [];
    % Iterate over components
    for i = I
        Pi = kron(P(:,:,i),V(:,:,i))/(v(i)+s-s*d-2);
        % Do not consider weights 求较大的项与最大项的距离
        [~,~,G_diff,N_diff,IW_diff] = ...
            GGIW_KLdiff(1,a(1,j),b(1,j),m(:,j),Pj,v(1,j),V(:,:,j),...
            1,a(1,i),b(1,i),m(:,i),Pi,v(1,i),V(:,:,i));
        

        if (G_diff <= U_G) && (N_diff <= U_N) && (IW_diff <= U_IW)
            L = [L i];
        end
    end
    
    if plotReduction
        for jj=L
            plotGGIWcomponent(a(jj),b(jj),m(:,jj),P(:,:,jj),v(jj),V(:,:,jj),'g',3)
        end
        title([num2str(length(L)) ' --- ' num2str(sum(w(1,L)))])
    end
    
    if length(L)>1
        % Save the merged weight
        wb(1,l) = sum(w(1,L));
        % Merge the components
        [wb(l),ab(l),bb(l),mb(:,l),Pb(:,:,l),vb(l),Vb(:,:,l),~] = ...
            GGIW_merge(w(L),a(L),b(L),m(:,L),P(:,:,L),v(L),V(:,:,L),length(L),model);
    else
        wb(1,l) = w(1,L);
        ab(l) = a(L);
        bb(l) = b(L);
        mb(:,l) = m(:,L);
        Pb(:,:,l) = P(:,:,L);
        vb(l) = v(L);
        Vb(:,:,l) = V(:,:,L);
    end
    
    if plotReduction
        plotGGIWcomponent(ab(l),bb(l),mb(:,l),Pb(:,:,l),vb(l),Vb(:,:,l),'k')
    end
    
    % Remove the merged components from the set of components
    I = setdiff(I,L);
end

% To get the same sum of weights
sumwb = sum(wb);
% Check if there are too many elements
Jb = min(l,Jmax);
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
vb = vb(1,i);
Vb = Vb(:,:,i);  
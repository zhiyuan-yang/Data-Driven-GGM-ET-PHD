function [ghat_k,Xhat_k,Phat_k,ExtHat_k] = GGIW_stateExtraction(w_k,a_k,b_k,m_k,P_k,nu_k,V_k,J_k)
% Function that extracts state estimates from the PHD filter.
%ghat_k 量测率
%Xhat_k 运动状态
%Phat_k 状态协方差
%ExtHat_k 扩展目标状态
    d = size(V_k,1);
    s = size(P_k,1);

    ghat_k = zeros(1,0);
    Xhat_k = zeros(s*d,0);
    Phat_k = zeros(s*d,s*d,0);
    ExtHat_k = zeros(d,d,0);
    counter = 0;
    for i = 1:J_k
        if w_k(1,i) > 0.5
            counter = counter+1;
            ghat_k = [ghat_k a_k(i)/b_k(i)];
            Xhat_k = [Xhat_k , m_k(:,i)];
            Phat_k(:,:,counter) = kron(P_k(:,:,i),V_k(:,:,i))/max(1,nu_k(1,i)+s-s*d-2);
            ExtHat_k(:,:,counter) = V_k(:,:,i)/max(1,(nu_k(1,i)-2*d-2));
        end
    end
end
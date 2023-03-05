function [ghat_k,Xhat_k,Phat_k,ExtHat_k] = ERHM_stateExtraction(w_k,a_k,b_k,m_k,P_k,J_k,model)
% Function that extracts state estimates from the PHD filter.
%ghat_k 量测率
%Xhat_k 运动状态
%Phat_k 状态协方差
%ExtHat_k 扩展目标状态
    d = model.RHM.d;
    x_dim = model.RHM.x_dim;
    
    ghat_k = zeros(1,0);
    Xhat_k = zeros(x_dim,0);
    Phat_k = zeros(x_dim,x_dim,0);
    ExtHat_k = zeros(d,d,0);
    counter = 0;
    for i = 1:J_k
        if w_k(1,i) > 0.5
            counter = counter+1;
            ghat_k = [ghat_k a_k(i)/b_k(i)];
            Xhat_k = [Xhat_k , m_k(1:x_dim,i)];
            Phat_k(:,:,counter) = P_k(1:x_dim,1:x_dim,i);
            L_k = [m_k(x_dim+1,i),              0;
                   m_k(x_dim+3,i),m_k(x_dim+2,i)];
            ExtHat_k(:,:,counter) = L_k*L_k';
        end
    end
end
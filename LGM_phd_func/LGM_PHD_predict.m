function [w_k,a_k,b_k,m_k,P_k,J_k] = LGM_PHD_predict(w_k1,alpha_k1,beta_k1,m_k1,P_k1,J_k1,model,ego)
%prediction process of phd filter
%input:   w_k1:       1*J_k1 vector
%         alpha_k1:   1*J_k1 vector
%         beta_k1     1*J_k1 vector
%         m_k1        dimx*J_k1 matrix
%         P_k1        dimx*dimx*J_k1 matrix
%         J_k1
%         model


%% temporal evoluation of the kinematic state
F_k1 = model.LGM.F_k;
Q_k1 = model.LGM.Q_k;
eta_k = model.LGM.eta_k;

%% 参数声明
% Probabilities of survival and detection
p_S_k = model.kin.pS;
% The spontaneous birth distributions
w_gam_k = model.birth.w;
a_gam_k = model.birth.a;
b_gam_k = model.birth.b;
m_gam_k = model.birth.m_LGM;
P_gam_k = model.birth.P_LGM;
J_gam_k = model.birth.J;

% The spawn distributions
J_beta_k = 0;% no spawn consider

%% 给需要预测的参数预分配内存
w_k = zeros(1,J_gam_k+J_beta_k*J_k1+J_k1);
a_k = zeros(1,J_gam_k+J_beta_k*J_k1+J_k1);
b_k = zeros(1,J_gam_k+J_beta_k*J_k1+J_k1);
m_k = zeros(size(m_k1,1),J_gam_k+J_beta_k*J_k1+J_k1);
P_k = zeros(size(P_k1,1),size(P_k1,2),J_gam_k+J_beta_k*J_k1+J_k1);

%% 新生目标
i = 0;                                                                     % 所有的预测项数 = 新生 + 存活
for j = 1:J_gam_k
    i = i+1;
    w_k(1,i) = w_gam_k(1,j);
    a_k(1,i) = a_gam_k;
    b_k(1,i) = b_gam_k;
    m_k(:,i) = m_gam_k(:,j) + [ego(1:2,:);zeros(4,1)];
    P_k(:,:,i) = P_gam_k;
end

%% 存活目标预测
for j = 1:J_k1
    i = i+1;
    w_k(1,i) = p_S_k*w_k1(1,j);
    a_k(1,i) = alpha_k1(1,j)/eta_k;                       
    b_k(1,i) = beta_k1(1,j)/eta_k;
    m_k(:,i) = F_k1*m_k1(:,j);                   
    P_k(:,:,i) = Q_k1 + F_k1*P_k1(:,:,j)*F_k1';       
    P_k(:,:,i) = 0.5*(P_k(:,:,i)+P_k(:,:,i)');                             % 使矩阵对称
end
J_k = i;

end
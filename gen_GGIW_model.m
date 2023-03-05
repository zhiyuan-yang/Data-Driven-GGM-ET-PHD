function model = gen_GGIW_model()
%% basic parameters
model.kin.pS = .99;                                                        % survival/death parameters, define p_S as a constant

%% Ŀ����״����
model.Ellipse.sigma_ex = 5;                                                    %������
model.Ellipse.sigma_ey = 15;
model.Ellipse.l_dim = 3;                                         % shape state dim for Ellipse


%% Ŀ���˶�ģ�� (GGIW model)
model.GGIW.x_dim = 4;                                                       % dimension of state vector,��nx=s*d��
model.GGIW.s = 2;
theta = 0.5;                                                                 % �������ʱ��
SIG = 0.1;                                                                 % ���ٶȱ�׼��
model.GGIW.d = 2;
model.GGIW.Ts = 0.1;                                                       % �������
model.GGIW.Id = eye(model.GGIW.d);
model.GGIW.F_k = [1,model.GGIW.Ts;0,1];
model.GGIW.Q_k = [model.GGIW.Ts^3/3,model.GGIW.Ts^2/2;model.GGIW.Ts^2/2,model.GGIW.Ts];         % ��������Э�������
model.GGIW.X_k = diag([25 25]);
model.GGIW.H_k = [1 0];                                                   % Measurement model


model.GGIW.v_min=model.GGIW.d+5;
we = 25;                                                                   % ����
model.GGIW.eta_k = we/(we-1);                                               % ��������
model.GGIW.tau = 5 ;                                                     % ʱ���ӳٳ���


%% birth parameters
model.birth.x_sym  = 20;                % birth positions 30*5 square around ego car
model.birth.y_sym  = 4;                                                       % number of births
model.birth.J = 4;                       % number of births
model.birth.w = [0.1 0.1 0.1 0.1];       % weights of birth components                                         % weights of birth components
model.birth.a = 15;
model.birth.b = 5;


model.birth.P_GGIW = diag([5 10]);
model.birth.v = model.GGIW.v_min;                                          % IW�ֲ����ɶ�v,must be >d+1 to compute likelihood, and predict nu, V
model.birth.V = diag([1 1]);   % birth covariances

model.max_card = 100;    %���Ŀ����
model.birth.card = zeros(1,model.max_card +1);
model.birth.card(1) = 0.90;                                 % Prob of 0 birth
model.birth.card(2) = 0.05;                                 % Prob of 1 birth
model.birth.card(3) = 0.025;                                % Prob of 2 birth
model.birth.card(4) = 0.025;                                % Prob of 3 birth
%% observation model prams
model.obs.z_dim = 2;                                                       % dimension of observation vector
model.obs.xrange = [-30 30];                                           % monitoring region dimensions in meters
model.obs.yrange = [-30 30];
model.obs.H = [1 0 0 0;
    0 1 0 0 ];                                               % observation matrix
model.obs.sigma_r = 2;                                                  % observation noise std in meters
model.obs.R = model.obs.sigma_r^2 * eye(model.obs.z_dim);                  % measurement noise covariance
model.obs.gamma = 50; % gamma Խ��RHM����Խ��
model.obs.pD = .99;                                                         % detection parameters, define p_D as a constant

%% clutter parameters
model.obs.lambda_c = 2;                                                    % poisson average rate of clutter (per scan) and identical for all sensors (this could be made sensor dependent in future releases)
model.obs.range_c  = [-20 20;-20 20];                              % uniform clutter on XY
model.obs.pdf_c    = 1/prod(model.obs.range_c(:,2)-model.obs.range_c(:,1));% uniform clutter density


%% OSPA parameters
model.ospa_cutoff = 100;
model.ospa_order  = 1;

%% Gate parameters, Distance probabilities
max_dist_prob = 0.80;
min_dist_prob = 0.40;
% Max and min distance for measurement set partitioning
model.obs.r =5;
model.max_distance = model.obs.r*chi2inv(max_dist_prob,model.obs.z_dim); % ���Ŷ�Ϊ0.8ʱ�ľ���
model.min_distance = model.obs.r*chi2inv(min_dist_prob,model.obs.z_dim); % ���Ŷ�Ϊ0.4ʱ�ľ���


%% For solving the equation Phat = kron(P,V) % Ϊ����˶�״̬Э�������
Pidx = 1:(model.GGIW.s * model.GGIW.d)^2;
Pidx = reshape(Pidx,model.GGIW.s * model.GGIW.d,model.GGIW.s * model.GGIW.d)'; % ��Pidx��дΪs*d��s*d�еľ���
Pidxvec = [];
for r=1:model.GGIW.s
    for c=1:model.GGIW.s
        Pidxsub = Pidx((1:model.GGIW.d)+model.GGIW.d*(r-1),(1:model.GGIW.d)+model.GGIW.d*(c-1))';
        Pidxvec = [Pidxvec;Pidxsub(:)];
    end
end
model.idxp = Pidxvec;

%% scenario1
model.ini_J = 2;                                                         % GGM ��������    
model.ini_w = [0.9,0.9];
model.ini_b = [5,5];                                               
model.ini_a = [20,20];                                                   % gamma�����ľ�ֵΪ�� a/b = 30/1 = beta_D1  
model.ini_m= [58 0 15 0;
              32 -4 15 0]';                                         % m = [x y vx vy l w]                                                                   
model.ini_P = diag([0.01 10]);
model.ini_P = cat(3,model.ini_P,model.ini_P);
model.ini_V = cat(3,eye(2),eye(2));
model.ini_nu = [7,7];

model.birth.m_GGIW = [-model.birth.x_sym  -model.birth.y_sym  15 0;
                      model.birth.x_sym  -model.birth.y_sym  15 0;
                     -model.birth.x_sym   model.birth.y_sym   15 0;
                      model.birth.x_sym   model.birth.y_sym   15 0]';  % birth means
model.lambda = 20;   %�����ⵥԪ��������������10ʱ�ָ�����

%% scenario2
% model.ini_J = 1;                                                           % GGM ��������    
% model.ini_w = 0.9;
% model.ini_b = 5;                                               
% model.ini_a = 40;                                                          % gamma�����ľ�ֵΪ�� a/b = 30/1 = beta_D1  
% model.ini_m= [45 -1.8 15 0]';                                              % m = [x y vx vy l w]                                                                   
% model.ini_P = diag([0.01 10]);
% model.ini_V = eye(2);
% model.ini_nu = 7;
% model.lambda = 11;   %�����ⵥԪ��������������10ʱ�ָ�����
% model.birth.m_GGIW = [model.birth.x_sym, -model.birth.y_sym, 0, 15; 
%                     -model.birth.x_sym, -model.birth.y_sym, 0, 15; 
%                     -model.birth.x_sym,  model.birth.y_sym  0 15;
%                     model.birth.x_sym, model.birth.y_sym, 0, 15]';  % birth means

end
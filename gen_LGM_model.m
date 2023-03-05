function model = gen_LGM_model()
%% basic parameters
model.kin.pS = .99;                                                        % survival/death parameters, define p_S as a constant
model.Ts = 0.1;
alpha = 0.0001; %扩展形状遗忘因子


%% LGM Model for RHM

model.LGM.Id = 2;
model.LGM.Ts = model.Ts;
model.LGM.x_dim = 6;    % x = [x y vx vy l w]^T;
model.LGM.kin_dim = 4;  % kinematic part
model.LGM.ex_dim = 2;   % extent part
model.LGM.F_k1 = [1 model.Ts;0 1];
model.LGM.F_ke = exp(-alpha * model.LGM.Ts)*eye(model.LGM.Id);
model.LGM.F_k = blkdiag(kron(model.LGM.F_k1,eye(model.LGM.Id)),model.LGM.F_ke);
model.LGM.Q_1 = kron([model.LGM.Ts^3/3,model.LGM.Ts^2/2;model.LGM.Ts^2/2,model.LGM.Ts] ...
                      ,blkdiag(0.01^2,0.01^2)); 
model.LGM.Q_e = (1-exp(-2*alpha*model.LGM.Ts))*eye(model.LGM.Id);
model.LGM.Q_k = blkdiag(model.LGM.Q_1,model.LGM.Q_e);
model.LGM.eta_k = 25/24;  %伽马分布参数遗忘因子

%% birth parameters
model.birth.x_sym  = 20;                % birth positions 30*5 square around ego car
model.birth.y_sym  = 4;
model.birth.J = 4;                       % number of births
model.birth.w = [0.1 0.1 0.1 0.1];       % weights of birth components
model.birth.a = 15;
model.birth.b = 5;


model.birth.P_LGM = diag([5 5 10 10 ones(1,2)]);

model.max_card = 100;    %最大目标数
model.birth.card = zeros(1,model.max_card +1);
model.birth.card(1) = 0.90;                                 % Prob of 0 birth
model.birth.card(2) = 0.05;                                 % Prob of 1 birth
model.birth.card(3) = 0.025;                                % Prob of 2 birth
model.birth.card(4) = 0.025;                                % Prob of 3 birth

%% observation model prams
model.obs.z_dim = 2;                                                       % dimension of observation vector
model.obs.xrange = [-30 30];                                           % monitoring region dimensions in meters
model.obs.yrange = [-30 30];
model.obs.H = [1 0 0 0 0 0;
               0 1 0 0 0 0];                                               % observation matrix
model.obs.sigma_r = 2;                                                     % observation noise std in meters
model.obs.R = model.obs.sigma_r^2 * eye(model.obs.z_dim);                  % measurement noise covariance
model.obs.gamma = 50; % gamma 越大，RHM优势越大
model.obs.pD = .99;                                                        % detection parameters, define p_D as a constant

%% clutter parameters
model.obs.lambda_c = 2;                                                    % poisson average rate of clutter (per scan) and identical for all sensors (this could be made sensor dependent in future releases)
model.obs.range_c  = [-20 20;-20 20];                                      % uniform clutter on XY
model.obs.pdf_c    = 1/prod(model.obs.range_c(:,2)-model.obs.range_c(:,1));% uniform clutter density

%% Plot parameters
model.plot.xrange = [-700 830];
model.plot.yrange = [-600 630];

%% OSPA parameters
model.ospa_cutoff = 100;
model.ospa_order  = 1;

%% Gate parameters, Distance probabilities
max_dist_prob = 0.80;
min_dist_prob = 0.40;
% Max and min distance for measurement set partitioning
model.obs.r = 5;
model.max_distance = model.obs.r*chi2inv(max_dist_prob,model.obs.z_dim); % 置信度为0.8时的距离
model.min_distance = model.obs.r*chi2inv(min_dist_prob,model.obs.z_dim); % 置信度为0.4时的距离


%% scenario1
model.ini_J = 2;                                                         % GGM 分量个数    
model.ini_w = [0.9,0.9];
model.ini_b = [5,5];                                               
model.ini_a = [20,20];                                                   % gamma函数的均值为： a/b = 30/1 = beta_D1  
model.ini_m= [58 0 15 0 4 2;
              32 -4 15 0 4 2]';                                         % m = [x y vx vy l w]                                                                   
model.ini_P = diag([0.01 0.01 10 10 10 10]);
model.ini_P = cat(3,model.ini_P,model.ini_P);

model.birth.m_LGM = [model.birth.x_sym  -model.birth.y_sym  15 0 4 2
                    -model.birth.x_sym  -model.birth.y_sym  15 0 4 2
                    -model.birth.x_sym   model.birth.y_sym  15 0 4 2;
                    model.birth.x_sym   model.birth.y_sym  15 0 4 2]';  % birth means
model.lambda = 15;   %当量测单元内量测数量超过10时分割量测


%% scenario2
% model.ini_J = 1;                                                           % GGM 分量个数    
% model.ini_w = 0.9;
% model.ini_b = 5;                                               
% model.ini_a = 40;                                                          % gamma函数的均值为： a/b = 30/1 = beta_D1  
% model.ini_m= [45 -1.8 15 0 4 2]';                                              % m = [x y vx vy l w]                                                                   
% model.ini_P = diag([0.01 0.01 10 10 10 10]);

% model.birth.m_LGM = [model.birth.x_sym  -model.birth.y_sym   0 15 4 2
%                     -model.birth.x_sym  -model.birth.y_sym  0 15  4 2
%                     -model.birth.x_sym   model.birth.y_sym  0 15  4 2;
%                     model.birth.x_sym   model.birth.y_sym  0 15  4 2]';  % birth means

% model.lambda = 11;   %当量测单元内量测数量超过10时分割量测

end


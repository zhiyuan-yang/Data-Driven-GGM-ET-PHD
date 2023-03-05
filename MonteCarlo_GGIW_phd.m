    %% setup
    clear;
    clc;
    addpath('GGIW_phd_func\');
    addpath('_common');
    model = gen_GGIW_model();                                               %load model
    filter = gen_GGIW_filters();                                            %load filter
    rng(2022); 
    montecarlo = 200;

    %scenario1 
%      load('data/scenario1/scenario1.mat');
%     %note: meas/truth/ego is under global coordinate system
%     meas = scenario1.meas;                                                 %1*k cell each cell is a 6*1 vector [x y z vx vy vz]
%     truth = scenario1.truth;                                               %1*k cell each cell is a struct with position 3*1 vctor and ges 3*1 vector
%     ego = scenario1.ego;                                                   %1*k cell each cell is a struct with position 3*1 vctor and ges 3*1 vector
    
%sccenario2
    load('data/scenario3/scenario3.mat'); 
    %note: meas/truth/ego is under global coordinate system
    meas = scenario3.meas;                                                 %1*k cell each cell is a 6*1 vector [x y z vx vy vz]
    truth = scenario3.truth;                                               %1*k cell each cell is a struct with position 3*1 vctor and ges 3*1 vector
    ego = scenario3.ego;  
    
    outcome = cell(montecarlo,1);
    %% output variables
    for i=1:montecarlo
        outcome{i}.output_GGIW_phd.state_X_hat = cell(length(meas),1);                          % inferred target states
        outcome{i}.output_GGIW_phd.state_E_hat = cell(length(meas),1);
        outcome{i}.output_GGIW_phd.state_n_hat = zeros(length(meas),1);                         % estimated number of targets
        outcome{i}.output_GGIW_phd.filter = filter;                                       % store filtering params
        outcome{i}.output_GGIW_phd.state_n_hat(1) = 2;
    end
    %% 初始化分布
    clear w a b m P nu V J                                                  % 用于GGIW_PHD初始化
    J = model.ini_J;                                                       % GGM 分量个数    
    w = model.ini_w;
    b = model.ini_b;                                               
    a = model.ini_a;                                                       % gamma函数的均值为： a/b = 30/1 = beta_D1  
    m = model.ini_m;                                                       % m = [x y vx vy l w]                                                                   
    P = model.ini_P; 
    nu = [7,7];                                                              % 自由度
    V = diag([1 1]);                                                       % 逆尺度矩阵
    V = cat(3,V,V);
    partition = cell(length(meas),1);                                      % 方差
  

    %% 滤波
    for i=1:montecarlo
        for k = 2:length(meas)
           fprintf('Current Monte Carlo Run: %d, Current time step: %d\n', i,k-1);
           %Predict
           [w,a,b,m,P,nu,V,J] = GGIW_PHD_predict(w,a,b,m,P,nu,V,J,model,ego{k}.position);

           %当前帧量测
           Zk = meas{k};
           Zkp = distance_partition_with_sub_partition(Zk,model,filter);
           partition{k} = Zkp;

           %updte
           [w,a,b,m,P,nu,V,J] = GGIW_PHD_update(w,a,b,m,P,nu,V,J,Zkp,model);

           %假设删减
           [w,a,b,m,P,nu,V,J] = GGIW_reduction(w,a,b,m,P,nu,V,J,model,filter);


            %状态提取
            [~,outcome{i}.output_GGIW_phd.state_X_hat{k},~,outcome{i}.output_GGIW_phd.state_E_hat{k}] = ...
                GGIW_stateExtraction(w,a,b,m,P,nu,V,J);
            outcome{i}.output_GGIW_phd.state_n_hat(k) = size(outcome{i}.output_GGIW_phd.state_X_hat{k},2);
        end
    end
    %% return
    % return filter computation time
    % output_GGIW_phd.all_time = toc(GGIW_phd_mark_start);
    ave_ospa = zeros(1,length(truth));
    for i=1:montecarlo
        % calculate and display average OSPA error
        error_t = zeros(1,length(truth));
        for k = 1:length(truth)
            error_t(k) = ospa_dist( get_comps(truth{k}.position,[1,2]), get_comps(outcome{i}.output_GGIW_phd.state_X_hat{k},[1 2]), ...
                model.ospa_cutoff, model.ospa_order);
            ave_ospa(k) = ave_ospa(k) + error_t(k);
        end
        outcome{i}.output_GGIW_phd.ospa = error_t;
    end
    ave_ospa = ave_ospa/montecarlo;
    t = 1:1:length(truth);
    figure(1);
    plot(t,ave_ospa);


%% function
function Xc= get_comps(X,c)
    if isempty(X)
        Xc= [];
    else
        Xc= X(c,:);
    end
end
    %% setup
    clear;
    clc;
    addpath('LGM_phd_func\');
    addpath('_common');
    model = gen_LGM_model();                                               
    filter = gen_LGM_filters();                                                                              
    load('data/GMModel.mat');

%     %% scenario1
    load('data/scenario1/scenario1.mat');
    %note: meas/truth/ego is under global coordinate system

    meas = scenario1.meas;                                                 %1*k cell each cell is a 6*1 vector [x y z vx vy vz]
    truth = scenario1.truth;                                               %1*k cell each cell is a struct with position 3*1 vctor and ges 3*1 vector
    ego = scenario1.ego;                                                   %1*k cell each cell is a struct with position 3*1 vctor and ges 3*1 vector

    %% scenario2
%     load('data/scenario3/scenario3.mat');
%     %note: meas/truth/ego is under global coordinate system
% 
%     meas = scenario3.meas;                                                 %1*k cell each cell is a 6*1 vector [x y z vx vy vz]
%     truth = scenario3.truth;                                               %1*k cell each cell is a struct with position 3*1 vctor and ges 3*1 vector
%     ego = scenario3.ego;                                                   %1*k cell each cell is a struct with position 3*1 vctor and ges 3*1 vector

    
    %% output variables1
    output_LGM_phd.state_X_hat = cell(length(meas),1);                     % inferred target states
    output_LGM_phd.state_E_hat = cell(length(meas),1);
    output_LGM_phd.state_n_hat = zeros(length(meas),1);                    % estimated number of targets
    output_LGM_phd.state_g_hat = cell(length(meas),1);                % estimated number of measurements of targets
    output_LGM_phd.filter = filter;                                        % store filtering params
    RHM_phd_mark_start = tic;                                              % start timing the filter
    
    %% ��ʼ���ֲ�
    clear w a b m P J                                                      % ����LGM_PHD��ʼ��
    J = model.ini_J;                                                       % GGM ��������    
    w = model.ini_w;
    b = model.ini_b;                                               
    a = model.ini_a;                                                       % gamma�����ľ�ֵΪ�� a/b = 30/1 = beta_D1  
    m = model.ini_m;                                                       % m = [x y vx vy l w]                                                                   
    P = model.ini_P; 
    output_LGM_phd.state_n_hat(1) = 2;
    partition = cell(length(meas),1);
       
    %% �˲�
    for k = 2:length(meas)
        fprintf('Current time step: %d\n', k-1);
        %Predict
        [w,a,b,m,P,J] = LGM_PHD_predict(w,a,b,m,P,J,model,ego{k}.position); 
        
        %��ǰ֡����
        Zk=meas{k};
        Zkp = distance_partition_with_sub_partition(Zk,model,filter);
        partition{k,1} = Zkp;

        %updte
        [w,a,b,m,P,J] = LGM_PHD_update(w,a,b,m,P,J,Zkp,model,filter,ego{k},GMModel);
        
        %����ɾ��
        [w,a,b,m,P,J] = LGM_reduction(w,a,b,m,P,J,filter);
        
        %״̬��ȡ
        [output_LGM_phd.state_g_hat{k},output_LGM_phd.state_X_hat{k},output_LGM_phd.state_P_hat{k},output_LGM_phd.state_E_hat{k}] = ...
                                                            LGM_stateExtraction(w,a,b,m,P,J,model);
        output_LGM_phd.state_n_hat(k) = size(output_LGM_phd.state_X_hat{k},2);
        
    end

    %% return
    % return filter computation time
    output_LGM_phd.all_time    = toc(RHM_phd_mark_start);

    % calculate and display average OSPA error
    error_t = zeros(1,length(truth));
    for k = 1:length(meas)
        error_t(k) = ospa_dist( get_comps(truth{k}.position,[1,2]), ....
            get_comps(output_LGM_phd.state_X_hat{k},[1 2]), model.ospa_cutoff, model.ospa_order);
    end
    output_LGM_phd.ospa = error_t;
    t = 1:1:length(meas);
    figure(1);
    plot(t,error_t);
      %% scenario1
    truth_n = 2*ones(1,length(truth));
    truth_n(1,20:80) = 1.5 * truth_n(1,20:80);
    

    %% scenario2
%     truth_n = 1*ones(1,length(truth));
%     truth_n(1,29:51) = 2 * truth_n(1,29:51);
    
    figure(2);
    hold on;
    plot(t,output_LGM_phd.state_n_hat);
    plot(t,truth_n);
    hold off;
%% function
function Xc= get_comps(X,c)
    if isempty(X)
        Xc= [];
    else
        Xc= X(c,:);
    end
end


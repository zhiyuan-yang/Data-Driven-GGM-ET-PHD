function filter_params = gen_ERHM_filters()
% Function that sets the parameters of several filters.

   
    %% ERHM_PHDfilter parameters
    filter_params.ERHM_phd.title              = 'ERHM-PHD ';
    filter_params.ERHM_phd.J_max              = 100;                       % maximum number of components
    filter_params.ERHM_phd.nb_comp_tg         = 4;                         % nb of components per target, used for pruning
    filter_params.ERHM_phd.T                  = 1e-2;                      % pruning threshold for components/hypotheses
    filter_params.ERHM_phd.U                  = 10;                           % mergng threshold for components/hypotheses
    filter_params.ERHM_phd.W_max              = 4;                         % maximum number of subsets
    filter_params.ERHM_phd.P_max              = 4;                         % maximum number of partitions
    filter_params.ERHM_phd.max_card           = 20;                          % maximum number of targets 
    filter_params.ERHM_phd.run_flag           = 'silence';                 %'disp' or 'silence' for on the fly output
    filter_params.N_low_FA           = 2;                           % % Number of measurements in a cell for it to be merged into FA cell
    filter_params.ERHM_phd.plot_partitions    = 0;                          % to plot the process of measurement
    filter_params.ERHM_phd.plot_process       = 0;                          % to plot the process of the predict and update
    filter_params.ERHM_phd.plot_reduction     = 0;                          % to plot the process of the reductio
    filter_params.ERHM_phd.alpha              = 1;                           % UT parameter
    filter_params.ERHM_phd.beta               = 2;
    filter_params.ERHM_phd.kappa              = 0;
    
  
end
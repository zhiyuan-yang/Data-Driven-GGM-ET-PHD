function filter_params = gen_GGIW_filters()
% Function that sets the parameters of several filters.

    %% GGIW_PHD(Kalman filter) filter parameters
    filter_params.GGIW_phd.title              = 'GGIW-PHD ';
    filter_params.GGIW_phd.J_max              = 100;                       % maximum number of components
    filter_params.GGIW_phd.nb_comp_tg         = 4;                         % nb of components per target, used for pruning
    filter_params.GGIW_phd.T                  = 1e-2;                      % pruning threshold for components/hypotheses
    filter_params.GGIW_phd.U_G                = 30;                         %¡¾"On the Reduction of Gaussian inverse Wishart Mixtures"¡¿GGIW-reduction
    filter_params.GGIW_phd.U_N                = 30;
    filter_params.GGIW_phd.U_IW               = 30;
    filter_params.GGIW_phd.W_max              = 4;                         % maximum number of subsets
    filter_params.GGIW_phd.P_max              = 4;                         % maximum number of partitions
    filter_params.GGIW_phd.max_card           = 20;                          % maximum number of targets 
    filter_params.GGIW_phd.run_flag           = 'silence';                 %'disp' or 'silence' for on the fly output
    filter_params.N_low_FA           = 2;                           % % Number of measurements in a cell for it to be merged into FA cell
    filter_params.GGIW_phd.plot_partitions    = 0;                          % to plot the process of measurement
    filter_params.GGIW_phd.plot_process       = 0;                          % to plot the process of the predict and update
    filter_params.GGIW_phd.plot_reduction     = 0;                          % to plot the process of the reduction

    
end
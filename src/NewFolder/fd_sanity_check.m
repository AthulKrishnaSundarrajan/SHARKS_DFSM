clc; clear; close all;

fun_name = 'two-link-robot';

switch fun_name

    case  'vanderpol'

        ex_output = @Vanderpol_output; % output function
        ex_plot = @Vanderpol_plot; % plot function
        
        
        n.ny = 2; % number of states
        n.nu = 1; % number of inputs
    
        % initialize
        f_og = cell(n.ny,1);
        
        f_og{1} = @(t,param,UYP) UYP(:,3);
        f_og{2} = @(t,param,UYP) -UYP(:,2)+UYP(:,3)-UYP(:,2).^2.*UYP(:,3)+UYP(:,1);
        
        % load optimal trajectory
        results = load('VDP_optimal_solution.mat');

    case 'two-link-robot'

        % number of controls, states, and parameters
        n.nu = 2; n.ny = 4; %n.np = 1;
        
        ex_output = @TwoLinkRobot_output; % output function
        ex_plot = @TwoLinkRobot_plot; % plot function
        
        % original derivative function
        f_og = cell(n.ny,1);
        f_og{1} = @(t,param,UYP) ((sin(UYP(:,5)).*(9/4*cos(UYP(:,5)).*UYP(:,3).^2+2*UYP(:,4).^2) + 4/3*(UYP(:,1)-UYP(:,2)) - 3/2*cos(UYP(:,5)).*UYP(:,2) )./ (31/36 + 9/4*sin(UYP(:,5)).^2));
        f_og{2} = @(t,param,UYP) (-( sin(UYP(:,5)).*(9/4*cos(UYP(:,5)).*UYP(:,4).^2+7/2*UYP(:,3).^2) - 7/3*UYP(:,2) + 3/2*cos(UYP(:,5)).*(UYP(:,1)-UYP(:,2)) )./ (31/36 + 9/4*sin(UYP(:,5)).^2) );
        f_og{3} = @(t,param,UYP) UYP(:,4) - UYP(:,3);
        f_og{4} = @(t,param,UYP) UYP(:,3);

        % load optimal trajectory
        results = load('TLR_optimal_results.mat');

end

% 'DFSM'

% define parameters
t0 = 0; tf = 5;

nt = 200; nsamples = 50;


% run simulation and get results
sim_details = run_simulation(t0,tf,nt,nsamples,fun_name);

% train-test split
split = [1,0];

% dfsm options
dfsm_options.ltype = 'LTI';
dfsm_options.ntype = 'GPR';
%dfsm_options.lsamples = 500;
dfsm_options.nsamples = 500;
dfsm_options.sampling_type = 'KM';
dfsm_options.train_test_split = split;
dfsm_options.scale_flag = ~true;

% construct dfsm model
dfsm = DFSM(sim_details,dfsm_options);

% extract dfsm elements
AB = dfsm.deriv.AB; % linear part

nonlin = dfsm.deriv.nonlin; % nonlinear part
error_ind = dfsm.deriv.error_ind; % error index

% convert the dfsm function into the form dtqp needs
f_dfsm = construct_function(AB,nonlin,error_ind,dfsm_options.ntype,n.ny,n.nu);

% extract results
T = results.T; U = results.U;Y = results.Y;

X = [U,Y];
param = [];

% evaluate derivative at the optimal trajectory
Df_of = DTQP_jacobian_real_central(f_og,X,T,param);
Df_dfsm = DTQP_jacobian_real_central(f_dfsm,X,T,param);

Df_diff = Df_of - Df_dfsm;

return
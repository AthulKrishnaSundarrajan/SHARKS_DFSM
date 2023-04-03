clc; clear; close all;

% number of controls, states, and parameters
n.nu = 2; n.ny = 4; %n.np = 1;

ex_output = @TwoLinkRobot_output; % output function
ex_plot = @TwoLinkRobot_plot; % plot function

dm_type = 'function';


switch dm_type

    case 'DFSM'

        % run simulations and construct DFSM for two-link-robot problem
        
        % define parameters
        t0 = 0; tf = 5;
        
        nt = 200; nsamples = 50;
        fun_name = 'two-link-robot';
        
        % run simulation and get results
        sim_details = run_simulation(t0,tf,nt,nsamples,fun_name);
        
        % train-test split
        split = [1,0];
        
        % dfsm options
        dfsm_options.ltype = '';
        dfsm_options.ntype = 'RBF';
        %dfsm_options.lsamples = 1000;
        dfsm_options.nsamples = 500;
        dfsm_options.sampling_type = 'KM';
        dfsm_options.train_test_split = split;
        dfsm_options.scale_flag = ~true;
        
        % construct dfsm model
        dfsm = DFSM(sim_details,dfsm_options);
         
        AB = dfsm.deriv.AB;
        
        nonlin = dfsm.deriv.nonlin;
        error_ind = dfsm.deriv.error_ind;

        f = construct_function(AB,nonlin,error_ind,dfsm_options.ntype,n.ny,n.nu);

        element.dynamics = []; % only needs to have the field to work
        dyn.f = f;
        setup.internalinfo.dyn = dyn;


    case 'function'

        % system dynamics
        f = cell(n.ny,1);
        f{1} = @(t,param,UYP) ((sin(UYP(:,5)).*(9/4*cos(UYP(:,5)).*UYP(:,3).^2+2*UYP(:,4).^2) + 4/3*(UYP(:,1)-UYP(:,2)) - 3/2*cos(UYP(:,5)).*UYP(:,2) )./ (31/36 + 9/4*sin(UYP(:,5)).^2));
        f{2} = @(t,param,UYP) (-( sin(UYP(:,5)).*(9/4*cos(UYP(:,5)).*UYP(:,4).^2+7/2*UYP(:,3).^2) - 7/3*UYP(:,2) + 3/2*cos(UYP(:,5)).*(UYP(:,1)-UYP(:,2)) )./ (31/36 + 9/4*sin(UYP(:,5)).^2) );
        f{3} = @(t,param,UYP) UYP(:,4) - UYP(:,3);
        f{4} = @(t,param,UYP) UYP(:,3);

        element.dynamics = []; % only needs to have the field to work
        dyn.f = f;
        setup.internalinfo.dyn = dyn;
end
%% tunable parameters
y0 = [0 0 0.5 0]; % initial states
yf = [0 0 0.5 0.522]; % final states

%% setup
% time horizon (scaled)
auxdata.t0 = 0; auxdata.tf = 1;


% Lagrange term
L(1).right = 1; L(1).left = 1; L(1).matrix = eye(2); % controls

% simple bounds
UB(1).right = 4; UB(1).matrix = y0; % initial states
LB(1).right = 4; LB(1).matrix = y0;
UB(2).right = 5; UB(2).matrix = yf; % final states
LB(2).right = 5; LB(2).matrix = yf;
UB(3).right = 1; UB(3).matrix = [1 1]; % controls
LB(3).right = 1; LB(3).matrix = [-1 -1];
%LB(4).right = 3; LB(4).matrix = 0.1; % parameters

% guess
Y0 = [y0;yf];
U0 = [[1 1];[-1 -1]];
%P0 = [[3.1];[3.1]];
%setup.guess.X = [U0,Y0,P0];
setup.guess.X = [U0,Y0];
% combine structures
setup.element = element;% setup.M = M; 
setup.UB = UB; setup.LB = LB; setup.L = L;
setup.t0 = auxdata.t0; setup.tf = auxdata.tf; setup.auxdata = auxdata; setup.n = n;

% solver opts

% opts
opts.general.displevel = 2;
opts.general.plotflag = 1;
opts.dt.defects = 'TR';
opts.dt.quadrature = 'CTR';
opts.dt.mesh = 'ED';
opts.dt.nt = 200; % number of nodes
opts.solver.display = 'iter'; % iterations
opts.solver.function = 'ipfmincon';
opts.method.form = 'nonlinearprogram';
opts.method.derivatives = 'real-central';
opts.solver.maxiters = 100;
opts.solver.tolerance = 1e-4;


%% solve
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

%% output
[O,sol] = ex_output(T,U,Y,P,F,in,opts);


%% plot
% disp("paused"); pause % for quasilinearization plots
ex_plot(T,U,Y,P,F,in,opts,sol)



return


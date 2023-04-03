clc; clear; close all;

dm_type = 'DFSM';

dfsm_type_cell = {'RBF','GPR'};

nstudies = length(dfsm_type_cell);

T_cell = cell(nstudies,1); Y_cell = cell(nstudies,1); U_cell = cell(nstudies,1);

ex_output = @Vanderpol_output; % output function
ex_plot = @Vanderpol_plot; % plot function


n.ny = 2; % number of states
n.nu = 1; % number of inputs


for istudy = 1:nstudies

    dfsm_type = dfsm_type_cell{istudy};
    
    switch dm_type
    
        case 'DFSM'
    
            % define parameters
            t0 = 0; tf = 5;
            
            nt = 200; nsamples = 50;
            fun_name = 'vanderpol';
            
            % run simulation and get results
            sim_details = run_simulation(t0,tf,nt,nsamples,fun_name);
            
            % train-test split
            split = [1,0];
            
            % dfsm options
            dfsm_options.ltype = 'LTI';
            dfsm_options.ntype = dfsm_type;
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
            f = construct_function(AB,nonlin,error_ind,dfsm_options.ntype,n.ny,n.nu);
            
            % combine
            element.dynamics = []; % only needs to have the field to work
            dyn.f = f;
            setup.internalinfo.dyn = dyn;
    
            
    
        case 'function'
            
            % initialize
            f = cell(n.ny,1);
    
            f{1} = @(t,param,UYP) UYP(:,3);
            f{2} = @(t,param,UYP) -UYP(:,2)+UYP(:,3)-UYP(:,2).^2.*UYP(:,3)+UYP(:,1);
    
            element.dynamics = [];
            dyn.f = f;
            setup.internalinfo.dyn = dyn;
    
    end
    
    
    
    L(1).left = 1; L(1).right = 1; L(1).matrix = 1;
    L(2).left = 2; L(2).right = 2; L(2).matrix = eye(2);
    setup.L = L;
    
    % simple bounds
    UB(1).right = 4; UB(1).matrix = [1;0]; % initial states
    LB(1).right = 4; LB(1).matrix = [1;0];
    UB(2).right = 1; UB(2).matrix = 1; % controls
    LB(2).right = 1; LB(2).matrix = -0.3;
    
    % guess
    Y0 = [[1,0];[1,0]];
    U0 = [[0];[0]];
    setup.guess.X = [U0,Y0];
    
    auxdata.t0 = 0; auxdata.tf = 5;
    
    setup.element = element; setup.UB = UB; setup.LB = LB;  setup.L = L;
    setup.t0 = auxdata.t0; setup.tf = auxdata.tf; setup.auxdata = auxdata; setup.n = n;
    
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
    
    T_cell{istudy} = T; U_cell{istudy} = U;Y_cell{istudy} = Y;
    
    %% output
    [O,sol] = ex_output(T,U,Y,P,F,in,opts);
    
    %% plot
    % disp("paused"); pause % for quasilinearization plots
    ex_plot(T,U,Y,P,F,in,opts,sol)

end

return


function sim_details = run_simulation(t0,tf,nt,nsamples,umin,umax,fun_name)

    % generate time
    time = linspace(t0,tf,nt);

    % generate samples
    state_names = {'x1','x2'};
    control_names = {'u'};
    
    % generate starting points
    nY0 = 2 + 2*rand(2,nsamples);
    
    % generate control samples
    u_samples = cell(nsamples,1);
    
    for i = 1:nsamples
        u_sample = lhsdesign_modified(nt,umin,umax);
    
        % create function
        u_pp = spline(time,u_sample);
        u_samples{i} = @(t) ppval(u_pp,t);
    
    
    end
    
    % initialize
    sim_details = cell(nsamples,1);
    
    % define solution options
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    
    % loop through and evaluate high fidelity function
    for isamples = 1:nsamples
    
        u_fun = u_samples{isamples};
        
        % run simulation
        [T,X] = ode45(@(t,y)og_deriv_function(t,y,u_fun,fun_name),[t0,tf],nY0(:,isamples),options);
    
        % evaluate control
        U = u_fun(T);
    
        % store
        sim_detail.time = T;
        sim_detail.states = X;
        sim_detail.controls = U;
        sim_detail.state_names = state_names;
        sim_detail.control_names = control_names;
        sim_detail.nstates = length(state_names);
        sim_detail.ncontrols = length(control_names);
        sim_detail.ninputs = length(state_names) + length(control_names);
        sim_detail.noutputs = length(state_names);
        sim_detail.fun_name = fun_name;
    
        
        % evaluate actual derivatives and store them
        dx_act = zeros(nt,2);
    
        for it = 1:length(T)
            % extract time and state value
            t = T(it);
            x = X(it,:);
            
            % evaluate original derivative function
            dx_act(it,:) = og_deriv_function(t,x,u_fun,fun_name);
        end
    
        % evaluate derivatives from polynomial approximation and store them
        x_pp = spline(T',X');
        dx_pp = fnder(x_pp,1);
        
        % evaluate polynomial approximate 
        dx_poly = ppval(dx_pp,T)';
    
        % store values
        sim_detail.state_derivatives = dx_poly;
        sim_detail.dx_act = dx_act;
    
        sim_details{isamples} = sim_detail;
    
    end
    
    % convert the cell of structure into a structure array
    sim_details = cell2struct_dfsm(sim_details,nsamples);

end
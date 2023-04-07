function sim_details = run_simulation(t0,tf,nt,nsamples,fun_name)

    % generate time
    time = linspace(t0,tf,nt);

    % generate samples
    switch fun_name

        case 'vanderpol'
            state_names = {'x1','x2'};
            control_names = {'u'};
            nx = 2;

            umin = -0.3; umax = 1;
            xmax = [5,5]; xmin = [-5,-5];
            x0 = [1,0];

            samp_type = 'LHC';

        case 'math2'
            state_names = {'x1','x2'};
            control_names = {'u'};
            nx = 2;

            umin = -5; umax = 5;

        case 'two-link-robot'
            state_names = {'x1','x2','x3','x4'};
            control_names = {'u1','u2'};
            nx = 4;

            umin = [-1,-1];umax = [1,1];
            xmax = ones(1,nx); xmin = -ones(1,nx);
            x0 = [0,0,0.5,0];
            samp_type = 'LHC';

        case 'transfer-min-fuel'

            state_names = {'r','t','vr','vt'};
            control_names = {'ur','ut'};
            nx = 4;

            umin = -[0.01,0.01]; umax = -umin;

            xmax = ones(1,nx); xmin = -ones(1,nx);
            x0 = [1,0,0,1];
            samp_type = 'freq-excite';




    end
    
    % generate starting points
    %nY0 = 2.*rand(nx,nsamples); 
    
    nY0 = lhsdesign_modified(nsamples,xmin,xmax);nY0 = nY0';

    nY0(:,end) = x0;
    
    % generate control samples
    u_samples = cell(nsamples,1);
    
    for i = 1:nsamples
        %u_sample = lhsdesign_modified(nt,umin,umax);
        u_sample = generate_usamples(samp_type,umin,umax,nt,time);
    
        % create function
        u_pp = spline(time,u_sample');
        u_samples{i} = @(t) ppval(u_pp,t);
    
    
    end
    
    % initialize
    sim_details = cell(nsamples,1);
    
    % define solution options
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);

    hf = figure;
    hf.Color = 'w';
    hold on;
    plot(nY0(1,:),nY0(2,:),'o')
    %xlim([xmin(1),xmax(1)]);ylim([xmin(2),xmax(2)]);
    
    % loop through and evaluate high fidelity function
    for isamples = 1:nsamples
    
        u_fun = u_samples{isamples};
        
        % run simulation
        [T,X] = ode45(@(t,y)og_deriv_function(t,y,u_fun,fun_name),[t0,tf],nY0(:,isamples),options);

        % plot states
        plot(X(:,1),X(:,2),'.','markersize',5)
    
        % evaluate control
        U = u_fun(T);
        
        if strcmpi(fun_name,'two-link-robot')
            U = U';
        end
    
        % store
        sim_detail.time = T;
        sim_detail.states = X;
        sim_detail.controls = U;
        sim_detail.state_names = state_names;
        sim_detail.control_names = control_names;
        sim_detail.output_names = [];
        sim_detail.outputs = [];
        sim_detail.nstates = length(state_names);
        sim_detail.ncontrols = length(control_names);
        sim_detail.ninputs = length(state_names) + length(control_names);
        sim_detail.nderiv = length(state_names);
        sim_detail.noutputs = [];
        sim_detail.fun_name = fun_name;
    
        
        % evaluate actual derivatives and store them
        dx_act = zeros(nt,nx);
    
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

function u_samples = generate_usamples(sampling_type,umin,umax,nt,time)
    
switch sampling_type

    case 'LHC'
        u_samples = lhsdesign_modified(nt,umin,umax);

    case 'freq-excite'
        nu = length(umin);
        u_samples = ones(nt,nu);
        perturb = 0.*rand(nt,nu);

        u_samples = u_samples.*(sin(exp(0.6.*time').*time')) + perturb;


end

end
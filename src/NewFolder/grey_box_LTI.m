clc; clear; close all;

rng(34534)

% file name

root_path = fileparts(which('INSTALL_DFSM'));
data_path = fullfile(root_path,'data');
fol_name_cell = 'FOWT_rated/16';



sim_path = fullfile(data_path,fol_name_cell);

% file names
prefix = 'IEA_w_TMD_';
suffix = '.outb';

% sim_files 
sim_dir  = dir(fullfile(sim_path,[prefix,'*',suffix]));

nsim = length(sim_dir);

if nsim <= 10
    numstring = '%01d';
else
    numstring = '%02d';
end


sim_files = cell(nsim,1);

% store the name of the files
for i = 1:nsim
    sim_files{i} = fullfile(sim_path,[prefix,num2str(i-1,numstring),suffix]);
end

% required states
%reqd_states = {'PtfmPitch','GenSpeed','YawBrTAxp'};

reqd_states = {'PtfmPitch','GenSpeed'};% ,'GenSpeed','YawBrTAxp'

% filtering arguments for states
filter_args.filt_states = [~true,~true,false];
filter_args.filt_states_tf = [1,1,1]; % 2,2 works

% required controls
reqd_controls = {'Wind1VelX','GenTq','BldPitch1','Wave1Elev'};
reqd_outputs = {'TwrBsMyt','GenSpeed','TwrBsFxt','YawBrTAxp','NcIMURAys','GenPwr'};%{'TwrBsFxt','TwrBsMyt','GenPwr','YawBrTAxp'};


% filtering arguments for controls
filter_args.filt_controls = [~true,false,false];
filter_args.filt_controls_tf = [1,0,0];

% filter flag
filter_flag = false;

%reqd_outputs = {};

% time
tmin = 00;
tmax = [];
add_dx2 = true;

% extract
sim_details = simulation_details(sim_files,reqd_states,reqd_controls,reqd_outputs,false,filter_flag,filter_args,add_dx2,tmin,tmax);


split = [0.3,0.7];

% dfsm options
dfsm_options.ltype = 'LTI';
dfsm_options.ntype = '';
dfsm_options.lsamples = nan;

dfsm_options.sampling_type = 'KM';
dfsm_options.train_test_split = split;
dfsm_options.scale_flag = 0;

nsamples = [30];
results_cell = cell(6,1);


dfsm_options.nsamples = 10;
% extract dfsm
%dfsm = DFSM(sim_details,dfsm_options);

% extract the linear model
% get the simulation data
[~,~,~,inputs,state_dx,outputs] = sample_data(sim_details(1:4),'KM',5);

AB = linsolve(inputs,state_dx)';
nc = length(reqd_controls);
ns = 2*length(reqd_states);
state_dx2 = state_dx(:,ns/2+1:end);

% initial parameter
par0 = AB(ns/2+1:end,:);
x0 = par0(:);
nx = length(x0);

%%

%
time = sim_details(1).time; dt = time(2) - time(1);
controls = inputs(:,1:nc);



%% testing out GA 
options_hybrid = optimoptions('fmincon','Display','iter','Algorithm','interior-point','MaxIterations',20,'MaxFunctionEvaluations',5000,'UseParallel',true);

% ga options
gaopt = optimoptions('ga','UseParallel',true,'Display','iter','MaxGenerations',nx*2);%,'HybridFcn',{@fmincon,options_hybrid});

% run hybrid optimization scheme
%x = ga(@(x) evaluate_loss(x,ns,nc,inputs,state_dx),nx,[],[],[],[],[],[],@(x)test_stability(x,ns,nc),[],gaopt) ; %@(x)findEIG(x,ns,nc)

% get the linear model

%opt_options.Advanced.StepTolerance = 1e-6;
opt_options.TolFun = 1e-6;
opt_options.TolX = 1e-6;
opt_options.MaxIter = 150;
opt_options.Algorithm = 'interior-point';

opt = greyestOptions;
opt.Display ='on';
opt.EnforceStability = true;
opt.SearchMethod ='fmincon';
opt.Focus = 'simulation';
opt.SearchOptions.MaxIterations = 15; 
opt.SearchOptions.FunctionTolerance = 1e-3;
%opt.OutputWeight = diag([1e-0,1e-0,1e-0,1e-0,1e-0,1e-0]);
opt.SearchOptions = opt_options;

% data


% initial system
% [A,B,C,D] = LTI_function(x,[],ns,nc);
% init_sys = ss(A,B,C,D);

% model
%m = idgrey('LTI_function',{'x',x},'c',{ns,nc});

%% organize data for multiple experiments
nt = length(time);

ne = length(outputs)/nt;

if ne == 1

    data = iddata(outputs,controls,dt);

elseif ne > 1

    data = cell(ne,1);

    
    for i = 1:ne

        [~,~,~,I,DX,Y] = sample_data(sim_details(nc),'KM',5);

        %Y(:,4:5) = Y(:,4:5)*1000;
        data{i} = iddata(Y,I(:,1:nc),dt);
    end

    data = merge(data{:});

end

%% N4SID

% initialize plot
hf = figure;
hf.Color = 'w';
hold on;
title('Poles','FontSize',10)
xlabel('Real','FontSize',10)
ylabel('Imag','FontSize',10)
axis equal; grid on
fplot(@(x) sqrt(1-x.^2),'color','k','LineWidth',1);
fplot(@(x) -sqrt(1-x.^2),'color','k','LineWidth',1);



nx_list = [4];
sys_cell = cell(length(nx_list),1);idx = 1;

func = 'ssest';
for nx_ = nx_list
    % run N4SID to get the system model and the starting point X0

    switch func

        case 'n4sid'

            % options
            options = n4sidOptions('InitialState','estimate','N4Weight','SSARX','Focus',...
                'simulation','EnforceStability',true,'Display','on');
        
            tic
            [sys,x0] = n4sid(data,'nx',nx_,'Ts',dt,options,'DisturbanceModel','none');
            t_train = toc;

        case 'ssest'

            % options
            options = ssestOptions('InitialState','estimate','N4Weight','SSARX','Focus',...
                'simulation','EnforceStability',true,'Display','on');

            tic
            [sys,x0] = ssest(data,1:16,'Ts',dt,options,'DisturbanceModel','none');
            t_train = toc;
    end



    
    % plot poles
    plot(real(eig(sys.A)),imag(eig(sys.A)),'r.','markersize',12);
    
    %% Test model
    % time and control inputs
    ind_test = 5;
    time = sim_details(ind_test).time;
    controls = sim_details(ind_test).controls;
    outputs = sim_details(ind_test).outputs;%outputs(:,4:5) = outputs(:,4:5)*1000;
    
    % starting point
    x0 = x0(:,1);
    u0 = controls(1,:);
    
    [nx,nu] = size(sys.B);
    ny = size(sys.C,1);
    
    X = zeros(nt+1,nx);
    Y = zeros(nt,ny);
    
    % extract state scape matrices
    A = sys.A;
    B = sys.B;
    C = sys.C;
    
    % X0
    X(1,:) = x0';

    tic;
    % loop through and evaluate response
    for i = 1:nt
        
        % input for the current time step
        u = controls(i,:)';
        
        % state value
        x = X(i,:)';
        
        % calculate state for next time step
        Xi = A*x + B*u;
        
        % calculate output
        Y(i,:) = C*Xi;

        % store state value for next iteration
        X(i+1,:) = Xi;       
    
    end
    t_eval = toc;
    % number of outputs
    ny = length(reqd_outputs);
    
    % loop through and plot outputs
    for i = 1:ny
    
        hf = figure;
        hf.Color = 'w';
        hold on;
    
        plot(time,outputs(:,i),'k-','LineWidth',1)
        plot(time,Y(:,i),'r-','LineWidth',1)
    
        xlim([time(1),time(end)])
        xlabel('Time [s]')
        legend('OpenFAST','N4SID')
        ylabel(reqd_outputs{i})
    
    end
    
    
    % loop throug and plots states
    for i = 1:nx
    
        hf = figure;
        hf.Color = 'w';
        hold on;
    
        plot(time,X(1:nt,i),'r-','LineWidth',1)
    
        xlim([time(1),time(end)])
        xlabel('Time [s]')
    
    end

    sys_cell{idx} = sys.A;
    idx = idx+1;

end
return
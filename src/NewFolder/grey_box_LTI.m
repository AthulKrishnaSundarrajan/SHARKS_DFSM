clc; clear; close all;

rng(34534)

% file name

root_path = fileparts(which('INSTALL_DFSM'));
data_path = fullfile(root_path,'data');
fol_name_cell = 'MHK_transition';



sim_path = fullfile(data_path,fol_name_cell);

% file names
prefix = 'RM1_';
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

reqd_states = {'PtfmSurge','PtfmSway','PtfmHeave','PtfmRoll','PtfmPitch','PtfmYaw','GenSpeed','YawBrTAxp'};% ,'GenSpeed','YawBrTAxp'

% filtering arguments for states
filter_args.filt_states = [~true,~true,false];
filter_args.filt_states_tf = [1,1,1]; % 2,2 works

% required controls
reqd_controls = {'RtVAvgxh','GenTq','BldPitch1','Wave1Elev'};

% filtering arguments for controls
filter_args.filt_controls = [~true,false,false];
filter_args.filt_controls_tf = [1,0,0];

% filter flag
filter_flag = false;

reqd_outputs = {};

% time
tmin = 00;
tmax = [];
add_dx2 = true;

% extract
sim_details = simulation_details(sim_files,reqd_states,reqd_controls,reqd_outputs,false,filter_flag,filter_args,add_dx2,tmin,tmax);

% extract the linear model
% get the simulation data
[~,~,~,inputs,state_dx,outputs] = sample_data(sim_details(1:2),'KM',5);

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
gaopt = optimoptions('ga','UseParallel',true,'Display','iter','MaxGenerations',nx*2,'HybridFcn',{@fmincon,options_hybrid});

% run hybrid optimization scheme
x = ga(@(x) evaluate_loss(x,ns,nc,inputs,state_dx),nx,[],[],[],[],[],[],@(x)test_stability(x,ns,nc),[],gaopt) ; %@(x)findEIG(x,ns,nc)

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
[A,B,C,D] = LTI_function(x,[],ns,nc);
init_sys = ss(A,B,C,D);

% model
m = idgrey('LTI_function',{'x',x},'c',{ns,nc});

% data
data = iddata(state_dx2,controls,dt);

% evaluate system
sys = greyest(data,m,opt);

% construct system
dfsm.deriv.AB = [sys.B,sys.A]';

ind_test = 3;
% time and control inputs
time = sim_details(ind_test).time;
controls = sim_details(ind_test).controls;


% test the dfsm

[dfsm_mf,X_cell,dx_cell,Y_cell] = test_dfsm(dfsm,sim_details(ind_test),ind_test,true,true);



return
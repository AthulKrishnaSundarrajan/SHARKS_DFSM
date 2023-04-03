clc; clear; close all;

% file name

root_path = fileparts(which('INSTALL_DFSM'));
data_path = fullfile(root_path,'data');
fol_name = '1p1_samples_belowrated';
sim_path = fullfile(data_path,fol_name);

% file names
prefix = 'lin_';
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
reqd_states = {'PtfmPitch','GenSpeed'};

% filtering arguments for states
filter_args.filt_states = [true,true];
filter_args.filt_states_tf = [2,2];

% required controls
reqd_controls = {'RtVAvgxh','GenTq','BldPitch1'};

% filtering arguments for controls
filter_args.filt_controls = [true,false,false];
filter_args.filt_controls_tf = [1,0,0];

% filter flag
filter_flag = true;

reqd_outputs = {'TwrBsFxt','TwrBsMxt'};

% time
tmin = 0;
tmax = [];
add_dx2 = true;

% extract
sim_details = simulation_details(sim_files,reqd_states,reqd_controls,reqd_outputs,filter_flag,filter_args,add_dx2,tmin,tmax);


split = [1,0];


% dfsm options
dfsm_options.ltype = 'LTI';
dfsm_options.ntype = 'RBF';
dfsm_options.lsamples = nan;
dfsm_options.nsamples = 500;
dfsm_options.sampling_type = 'KM';
dfsm_options.train_test_split = split;
dfsm_options.scale_flag = 0;

% extract dfsm
dfsm = DFSM(sim_details,dfsm_options);

% linear part
AB = dfsm.deriv.AB;

% nonlinear part
nonlin = dfsm.deriv.nonlin;
error_ind = dfsm.deriv.error_ind;

% number of states and controls
nx = sim_details(1).nstates;
nu = sim_details(1).ncontrols;

% change the form into something DTQP can use
f = construct_function(AB,nonlin,error_ind,dfsm_options.ntype,nx,nu);


element.dynamics = []; % only needs to have the field to work
dyn.f = f;
setup.internalinfo.dyn = dyn;

% opts
opts.general.displevel = 2;
opts.general.plotflag = 1;
opts.dt.defects = 'TR';
opts.dt.quadrature = 'CTR';
opts.dt.mesh = 'ED';
opts.dt.nt = 1000; % number of nodes
opts.solver.display = 'iter'; % iterations
opts.solver.function = 'ipfmincon';
opts.method.form = 'nonlinearprogram';
opts.method.derivatives = 'real-central';
opts.solver.maxiters = 250;
opts.solver.tolerance = 1e-5;

n.ny = nx; % number of states
n.nu = nu; % number of inputs

r1 = 1e-8; 
w1 = 1e-7; w2 = 1e-8;

% objective function
lx = 1;
% L(lx).left = 1;
% L(lx).right = 1;
% L(lx).matrix = diag([0,r1,0]);
% lx = lx+1;
% 
% L(lx).left = 1;
% L(lx).right = 2;
% Lmat = zeros(nu,nx); Lmat(2,2) = 1;
% L(lx).matrix = -Lmat;
% lx = lx+1;
%

% objective function

% min (19.8-tg)^2
lx = 1;
% L(lx).left = 0;
% L(lx).right = 0;
% L(lx).matrix = 19.8^2;
% 
% lx = lx+1;
L(lx).left = 1;
L(lx).right = 1;
L(lx).matrix = diag([0,1e-8,1e-0]);

% lx = lx+1;
% L(lx).left = 0;
% L(lx).right = 1;
% L(lx).matrix = [0,-2*19.8,0];

lx = lx+1;
%
L(lx).left = 1;
L(lx).right = 2;
Lmat = zeros(nu,nx); Lmat(2,2) = 1;
L(lx).matrix = -Lmat;
lx = lx+1;


% extract wind speed and construct interpolating function
controls = sim_details(2).controls;
states = sim_details(2).states;
time = sim_details(2).time;


x0 = states(1,:);
t0 = time(1); tf = time(end);

wind_speed = controls(:,1);

wind_pp = spline(time,wind_speed);
w_fun = @(t) ppval(wind_pp,t);


%% Control Bounds
ix = 1;

% controls upper bound
UB(ix).right = 1;
UB(ix).matrix = {@(t) w_fun(t);
    19.8;
    22};

% control lower bound
LB(ix).right = 1;
LB(ix).matrix = {@(t) w_fun(t);
    5;
    0};

ix = ix+1;

% states upper bound
UB(ix).right = 2;
UB(ix).matrix = [10,8,inf,inf];

% states lower bound
LB(ix).right = 2;
LB(ix).matrix = [0,2,-inf,-inf];

% states starting point
% ix = ix+1;
% UB(ix).right = 4;
% UB(ix).matrix = x0;
% 
% LB(ix).right = 4;
% LB(ix).matrix = x0;

% time span
auxdata.t0 = t0; auxdata.tf = tf;

% combine
setup.element = element; setup.UB = UB; setup.LB = LB;  setup.L = L;
setup.t0 = auxdata.t0; setup.tf = auxdata.tf; setup.n = n;

%% solve
[T,U,X,P,F,in,opts] = DTQP_solve(setup,opts);

close all;
% plot inputs
hf = figure;
hf.Color = 'w';
sgtitle('Controls')


for idx = 1:nu
    subplot(nu,1,idx)
    plot(T,U(:,idx),"LineWidth",1)
    xlabel('Time [s]')
    ylabel(sim_details(1).control_names{idx})
end

% plot
hf = figure;
hf.Color = 'w';
hold on;
sgtitle('State')

state_bounds = {[0,10],[0,8]};

for idx = 1:nx-2
    subplot(nx-2,1,idx)
    hold on;
    plot(T,X(:,idx),"LineWidth",1)
    xlabel('Time [s]')
    ylim(state_bounds{idx})
    ylabel(sim_details(1).state_names{idx})
end

return




clc; clear; close all;

% file name

root_path = fileparts(which('INSTALL_DFSM'));
data_path = fullfile(root_path,'data');
fol_name = 'DFSM_rated_10';
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
dfsm_options.ntype = 'GPR';
%dfsm_options.lsamples = 500;
dfsm_options.nsamples = 500;
dfsm_options.sampling_type = 'KM';
dfsm_options.train_test_split = split;
dfsm_options.scale_flag = 0;

dfsm = DFSM(sim_details,dfsm_options);

controls_mean = mean(sim_details(1).controls,1);
time = sim_details(1).time;

nt = length(time);

controls = ones(nt,length(controls_mean)).*controls_mean;

% t_ind = time >= 100 & time <= 200;
% controls(t_ind,3) = controls(t_ind,3) + 2;

t_ind = time >= 100;
controls(t_ind,3) = controls(t_ind,3)-1;

u_pp = spline(time,controls');
u_fun = @(t) ppval(u_pp,t);

controls1 = ones(nt,length(controls_mean)).*controls_mean;

u_pp1 = spline(time,controls1');
u_fun1 = @(t) ppval(u_pp1,t);

controls2 = ones(nt,length(controls_mean)).*controls_mean;
controls2(t_ind,3) = controls2(t_ind,3)+1.5;
u_pp2 = spline(time,controls2');
u_fun2 = @(t) ppval(u_pp2,t);
%------------------------------------------------------------------

x0 = sim_details(1).states(1,:);

[T,X] = ode45(@(t,x) ode_dfsm(t,x,u_fun,dfsm),[0,600],x0,odeset('RelTol',1e-6,'AbsTol',1e-6));
[T1,X1] = ode45(@(t,x) ode_dfsm(t,x,u_fun1,dfsm),[0,600],x0,odeset('RelTol',1e-6,'AbsTol',1e-6));
[T2,X2] = ode45(@(t,x) ode_dfsm(t,x,u_fun2,dfsm),[0,600],x0,odeset('RelTol',1e-6,'AbsTol',1e-6));
close all;
% plot inputs
hf = figure;
hf.Color = 'w';
sgtitle('Controls')

nc = size(controls,2);

for idx = 1:nc
    subplot(nc,1,idx)
    hold on
    plot(time,controls(:,idx),"LineWidth",1)
    plot(time,controls1(:,idx),"LineWidth",1)
    plot(time,controls2(:,idx),"LineWidth",1)
    xlabel('Time [s]')
    ylabel(sim_details(1).control_names{idx})
    legend('mean-1','mean','mean+1.5')
end

% plot states
hf = figure;
hf.Color = 'w';
hold on;
sgtitle('State')
nx = size(X,2);

for idx = 1:nx-2
    subplot(nx-2,1,idx)
    hold on;
    plot(T,X(:,idx),"LineWidth",1)
    plot(T1,X1(:,idx),"LineWidth",1)
    plot(T2,X2(:,idx),"LineWidth",1)
    legend('mean-1','mean','mean+1.5')
    xlabel('Time [s]')
    ylabel(sim_details(1).state_names{idx})
end

return
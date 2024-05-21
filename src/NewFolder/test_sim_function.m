clc; clear; close all;

%---------------------------------------
% test_sim_function.m
% A script to compare accuracy and speed between lsim and ode45 for
% simulating linear system
%---------------------------------------

% path to the data
root_path = fileparts(which('INSTALL_DFSM'));
data_path = fullfile(root_path,'data'); % data folder
fol_name_cell = 'DFSM_1p6_3';
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
reqd_outputs = {'TwrBsMyt','TwrBsFxt','NcIMURAys','GenPwr'};%{'TwrBsFxt','TwrBsMyt','GenPwr','YawBrTAxp'};


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


split = [0.5,0.5];

% dfsm options
dfsm_options.ltype = 'LTI';
dfsm_options.ntype = '';
dfsm_options.lsamples = nan;

dfsm_options.sampling_type = 'KM';
dfsm_options.train_test_split = split;
dfsm_options.scale_flag = 0;

nsamples = [5];
results_cell = cell(6,1);


dfsm_options.nsamples = 10;

% number of states and controls
ns = length(reqd_states)*2;
nc = length(reqd_controls);

% construct DFSM
%dfsm = DFSM(sim_details,dfsm_options);

[~,~,~,inputs,state_dx,outputs] = sample_data(sim_details(1:5),'KM',5);

% value to scale states by
states_max = max(abs(inputs(:,nc+1:end)));
states_max = states_max(1:ns/2);

% extract the linear systems
AB = linsolve(inputs,state_dx)';

% predicted derivatives
dx_pred = inputs*AB';

% error
error_dx = state_dx(:,ns/2+1:end) - dx_pred(:,ns/2+1:end);

% loss
loss_dx = 1/(length(error_dx))*trace(error_dx'*error_dx);

% extract A and B matrices
A = AB(:,nc+1:end);
B = AB(:,1:nc);

n_exp = 5;
t_ode45 = zeros(n_exp,1);
t_lsim = zeros(n_exp,1);
loss_exp = zeros(n_exp,1);


for i_exp = 1:n_exp

    % extract test indices
    [~,~,~,inputs,state_dx,outputs] = sample_data(sim_details(i_exp),'KM',5);
    
    % extract controls and states
    time = sim_details(i_exp).time;
    controls = inputs(:,1:nc);
    states = inputs(:,nc+1:end);
    
    
    tf = max(time);
    t_test_start = tf - 100;

    t_ind = (time >= t_test_start) & (time <= tf);

    time = time(t_ind);
    controls = controls(t_ind,:);
    states = states(t_ind,:);
    states_act = states(:,1:ns/2);
    
    % construct interpolating function for controls
    u_pp = interp1(time,controls,'pchip','pp');
    u_fun = @(t) ppval(u_pp,t);
    
    % initial states
    X0 = states(1,:);
    
    % run ode45
    tic;
    [T,X] = ode45(@(t,x)odefun(t,x,A,B,u_fun),time,X0);
    t_ode45(i_exp) = toc;
    
    % LSIM
    C = [eye(ns/2),zeros(ns/2)];
    D = zeros(ns/2,nc);
    
    % construct state space
    sys = ss(A,B,C,D);
    
    % simulate using lsim
    tic
    X_lsim = lsim(sys,controls,time,X0,'zoh');
    t_lsim(i_exp) = toc;
    
    
    % plot
    hf = figure;
    hf.Color = 'w';
    
    for i = 1:ns/2
    
        subplot(ns/2,1,i)
        hold on;
        plot(time,states_act(:,i),'LineWidth',1)
        plot(time,X(:,i),'LineWidth',1)
        plot(time,X_lsim(:,i),'LineWidth',1)
        xlim([time(1),time(end)]);
        ylabel(reqd_states{i});
      
        legend('OpenFAST','ode45','lsim')
    
    end
    
    xlabel('Time')
    
    % calculate error between actual and predicted state simulations
    SA = (states_act./states_max);
    SP = (X_lsim./states_max);
    error =  SA - SP; 
    
    % evaluate loss
    loss_exp(i_exp) = 1/(length(time))*(trace(error'*error));


end

t_lsim

t_ode45

loss_sim = mean(loss_exp);

return




function dx = odefun(t,x,A,B,u_fun)
    
    % controls
    u = u_fun(t);
    
    % evaluate dx
    dx = A*x + B*u;

end
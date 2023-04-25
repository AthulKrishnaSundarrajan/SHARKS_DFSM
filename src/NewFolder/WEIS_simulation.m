clc; clear; close all;

% file name

root_path = fileparts(which('INSTALL_DFSM'));
data_path = fullfile(root_path,'data');
fol_name = 'DFSM_transition_10';
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

reqd_outputs = {'TwrBsFxt'};

% time
tmin = 0;
tmax = [];
add_dx2 = false;

% extract
sim_details = simulation_details(sim_files,reqd_states,reqd_controls,reqd_outputs,filter_flag,filter_args,add_dx2,tmin,tmax);

split = [1,0];

% dfsm options
dfsm_options.ltype = 'LTI';
dfsm_options.ntype = 'GPR';
dfsm_options.lsamples = nan;
dfsm_options.nsamples = 500;
dfsm_options.sampling_type = 'KM';
dfsm_options.train_test_split = split;
dfsm_options.scale_flag = 0;

% extract dfsm
dfsm = DFSM(sim_details(1:2),dfsm_options);

[dfsm_mf,X_cell_mf,~] = test_dfsm(dfsm,sim_details(1),1,~false,~false);


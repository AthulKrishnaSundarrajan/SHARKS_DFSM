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

% time
tmin = 0;
tmax = [];

% extract
sim_details = simulation_details(sim_files,reqd_states,reqd_controls,filter_flag,filter_args,tmin,tmax);


return




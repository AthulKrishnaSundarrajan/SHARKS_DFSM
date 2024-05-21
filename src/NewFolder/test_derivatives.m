clc;clear;%close all;

rng(34534)

% file name

root_path = fileparts(which('INSTALL_DFSM'));
data_path = fullfile(root_path,'data');
fol_name_cell = 'transition';



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

reqd_states = {'PtfmPitch','TTDspFA','GenSpeed'};% ,'GenSpeed','YawBrTAxp'

% filtering arguments for states
filter_args.filt_states = [true,true,~false];
filter_args.filt_states_tf = 0.1*[1,1,1]; % 2,2 works

% required controls
reqd_controls = {'RtVAvgxh','GenTq','BldPitch1','Wave1Elev'};
nc = length(reqd_controls);
reqd_outputs = {};%{'TwrBsMyt','TwrBsFxt','YawBrTAxp','NcIMURAys','GenPwr'};%{'TwrBsFxt','TwrBsMyt','GenPwr','YawBrTAxp'};


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
time = sim_details(1).time;
[~,~,~,inputs,state_dx,outputs] = sample_data(sim_details(1:4),'KM',5);
[~,~,~,inputs1,state_dx1,outputs1] = sample_data(sim_details(1),'KM',5);

AB = linsolve(inputs,state_dx);

AB = AB';

A = AB(:,nc+1:end);

states = inputs1(:,nc+1:end);

hf = figure;
hf.Color = 'w';

tf = [0,0.1];

for tf_ = tf

state = perform_filter(states(:,1),time,tf_);

dt = time(2)-time(1);
st_pp = spline(time,state);
st_pp_dx = fnder(st_pp,1);
st_pp_dx2 = fnder(st_pp,2);
state_dx = ppval(st_pp_dx,time);
state_dx2 = ppval(st_pp_dx2,time);



subplot(3,1,1)
hold on;
plot(time,state,'-')

subplot(3,1,2)
hold on;
plot(time,state_dx,'-')

subplot(3,1,3)
hold on;
plot(time,state_dx2,'-')

end
return

function signal = perform_filter(signal,time,t_f)

if t_f == 0
    signal = signal;

else


    dt = time(2) - time(1);

    %
    nb = floor(t_f/dt);
    
    %
    b = ones(nb,1)/nb;
    
    %
    signal = filtfilt(b,1,signal);

end

end

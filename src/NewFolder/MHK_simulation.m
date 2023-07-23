clc; clear; close all;

% file name

root_path = fileparts(which('INSTALL_DFSM'));
data_path = fullfile(root_path,'data');
fol_name_cell = {'DFSM_MHK2'};

results_cell = cell(3,1);
dfsm_cell = cell(3,1);

for ifol = 1:length(fol_name_cell)

    sim_path = fullfile(data_path,fol_name_cell{ifol});
    
    % file names
    prefix = 'MHK_RM1_Floating_';
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
    reqd_states = {'PtfmSurge','PtfmSway','PtfmHeave','PtfmRoll','PtfmPitch','PtfmYaw'};
    
    % filtering arguments for states
    filter_args.filt_states = [~true,~true];
    filter_args.filt_states_tf = [1,1]; % 2,2 works
    
    % required controls
    reqd_controls = {'Wave1Elev'};
    
    % filtering arguments for controls
    filter_args.filt_controls = [~true,false,false];
    filter_args.filt_controls_tf = [1,0,0];
    
    % filter flag
    filter_flag = true;
    
    reqd_outputs = {};
    
    % time
    tmin = 0;
    tmax = [];
    add_dx2 = false;
    
    % extract
    sim_details = simulation_details(sim_files,reqd_states,reqd_controls,reqd_outputs,false,filter_flag,filter_args,add_dx2,tmin,tmax);
    
    split = [0.5,0.5];
    
    % dfsm options
    dfsm_options.ltype = '';
    dfsm_options.ntype = 'RBF';
    dfsm_options.lsamples = nan;
    dfsm_options.nsamples = 50;
    dfsm_options.sampling_type = 'KM';
    dfsm_options.train_test_split = split;
    dfsm_options.scale_flag = 0;
    
    % extract dfsm
    dfsm = DFSM(sim_details,dfsm_options);

    ind_test = 2;

    time = sim_details(ind_test).time;
    controls = sim_details(ind_test).controls;
    
    [dfsm,X_cell,~,Y_cell] = test_dfsm(dfsm,sim_details(ind_test),1,~false,~false);

    results_cell = {time,controls,X_cell,Y_cell};
    dfsm_cell{ifol} = dfsm;


end


input_train = sim_details(1).controls;

mat_name = 'DFSM_MHK_validation_EAB.mat';

save(mat_name,'results_cell','input_train')

return

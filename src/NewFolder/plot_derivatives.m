clc;clear; close all;

% set seed
rng('default')

% define parameters
t0 = 0; tf = 3;

nt = 200; nsamples = 100;
fun_name = 'two-link-robot';

% run simulation and get results
sim_details = run_simulation(t0,tf,nt,nsamples,fun_name);

split = [1,0];


%% train multifidelity dfsm
% dfsm options
dfsm_options.ltype = 'LTI';
dfsm_options.ntype = 'GPR';
dfsm_options.nsamples = 500;
dfsm_options.sampling_type = 'KM';
dfsm_options.train_test_split = split;
dfsm_options.scale_flag = ~true;


dfsm_mf = DFSM(sim_details,dfsm_options);

%% train a single nonlinear dfsm
% dfsm options
dfsm_options.ltype = '';
dfsm_options.ntype = 'GPR';
dfsm_options.nsamples = 500;
dfsm_options.sampling_type = 'KM';
dfsm_options.train_test_split = split;
dfsm_options.scale_flag = ~true;

% construct DFSM with the nonlinear parrt only
dfsm_nl = DFSM(sim_details,dfsm_options);

% initialize storage cell
dx_error_nl = cell(nsamples,1);
dx_error_mf = cell(nsamples,1);


% loop through and evaluate the error in the state derivatives predicted
ind = 1:nsamples;

[dfsm_mf,X_cell_mf,dx_cell_mf] = test_dfsm(dfsm_mf,sim_details,ind,false);

[dfsm_nl,X_cell_nl,dx_cell_nl] = test_dfsm(dfsm_nl,sim_details,ind,false);




% save results
mat_name = 'TLR_dx_error.mat';

save(mat_name,"dx_cell_nl","dx_cell_mf",'X_cell_nl','X_cell_mf');



%  
return



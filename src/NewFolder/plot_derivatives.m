clc;clear; close all;

% set seed
rng(4357)
fun_name = 'two-link-robot2';

% define parameters
t0 = 0; tf = 0.5;
fname = mfilename('fullpath');


fname = which(fname);

nt = 90; nsamples = 100;

% run simulation and get results
fac = 1;
x0 = zeros(1,4);
u0 = [0,0];
sim_details = run_simulation(t0,tf,nt,nsamples,fun_name,fac,x0,u0);

split = [0.8,0.2];

test_simulations_ind = floor(nsamples*0.8+1):nsamples;
test_simulations = sim_details(test_simulations_ind);


num_samples = [500];
construct_time = cell(length(num_samples),2);

for i = 1:length(num_samples)

%% train a single nonlinear dfsm
% dfsm options
dfsm_options.ltype = '';
dfsm_options.ntype = 'RBF';
dfsm_options.nsamples = num_samples(i);
dfsm_options.sampling_type = 'KM';
dfsm_options.train_test_split = split;
dfsm_options.scale_flag = ~true;


% construct DFSM with the nonlinear parrt only
dfsm_nl = DFSM(sim_details,dfsm_options);

%% train multifidelity dfsm
% dfsm options
dfsm_options.ltype = 'LTI';
dfsm_options.ntype = 'RBF';
dfsm_options.nsamples = num_samples(i);
dfsm_options.sampling_type = 'KM';
dfsm_options.train_test_split = split;
dfsm_options.scale_flag = ~true;
dfsm_options.use_samples.input_sampled = dfsm_nl.input_sampled;
dfsm_options.use_samples.dx_sampled = dfsm_nl.dx_sampled;
dfsm_options.use_samples.output_sampled = dfsm_nl.output_sampled;


dfsm_mf = DFSM(sim_details,dfsm_options);


%% train a single linear dfsm
% dfsm options
dfsm_options.ltype = 'LTI';
dfsm_options.ntype = '';
dfsm_options.nsamples = num_samples(i);
dfsm_options.sampling_type = 'KM';
dfsm_options.train_test_split = split;
dfsm_options.scale_flag = ~true;
dfsm_options.use_samples.input_sampled = dfsm_mf.input_sampled;
dfsm_options.use_samples.dx_sampled = dfsm_mf.dx_sampled;
dfsm_options.use_samples.output_sampled = dfsm_mf.output_sampled;

% construct DFSM with the nonlinear parrt only
dfsm_lin = DFSM(sim_details,dfsm_options);

construct_time{i,1} = sum(dfsm_nl.nonlin_construct);
construct_time{i,2} = sum(dfsm_mf.nonlin_construct) + dfsm_mf.lin_construct;


dfsm_options = rmfield(dfsm_options,'use_samples');
end
% loop through and evaluate the error in the state derivatives predicted
ind = randsample(1:0.2*nsamples,1);

time = sim_details(1).time;
%states = sim_details(ind).states;
controls = sim_details(1).controls;
%state_derivatives = sim_detials(ind).state_derivatives;

[dfsm_mf,X_cell_mf,~] = test_dfsm(dfsm_mf,test_simulations,1:0.2*nsamples,false,~false);

[dfsm_nl,X_cell_nl,~] = test_dfsm(dfsm_nl,test_simulations,1:0.2*nsamples,false,~false);

[dfsm_lin,X_cell_lin,~] = test_dfsm(dfsm_lin,test_simulations,1:0.2*nsamples,false,~false);

[dfsm_mf,~,dx_cell_mf] = test_dfsm(dfsm_mf,test_simulations,1:0.2*nsamples,false,false);

[dfsm_nl,~,dx_cell_nl] = test_dfsm(dfsm_nl,test_simulations,1:0.2*nsamples,false,false);

[dfsm_lin,~,dx_cell_lin] = test_dfsm(dfsm_lin,test_simulations,1:0.2*nsamples,false,false);



% save resul
mat_name = 'TLR_simulations.mat';

save(mat_name,"dx_cell_nl","dx_cell_mf",'dx_cell_lin','X_cell_nl','X_cell_mf','X_cell_lin','time','controls',"ind");



%  
return



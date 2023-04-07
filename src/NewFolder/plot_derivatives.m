clc;clear; close all;

% set seed
rng('default')

% define parameters
t0 = 0; tf = 1;

nt = 200; nsamples = 100;
fun_name = 'two-link-robot';

% run simulation and get results
sim_details = run_simulation(t0,tf,nt,nsamples,fun_name);

split = [0.9,0.1];

% dfsm options
dfsm_options.ltype = 'LTI';
dfsm_options.ntype = 'GPR';
dfsm_options.nsamples = 500;
dfsm_options.sampling_type = 'KM';
dfsm_options.train_test_split = split;
dfsm_options.scale_flag = ~true;


dfsm = DFSM(sim_details,dfsm_options);

X_sampled = vertcat(sim_details(:).states);
X0_sampled = zeros(nsamples,sim_details(1).nstates);

for i = 1:nsamples
    X0_sampled(i,:) = sim_details(i).states(1,:);
end

input_sampled = dfsm.input_sampled;
dx_sampled = dfsm.dx_sampled;

mat_name = 'TLR_state_samples.mat';

save(mat_name,'X_sampled','X0_sampled','input_sampled','dx_sampled');

%  
return



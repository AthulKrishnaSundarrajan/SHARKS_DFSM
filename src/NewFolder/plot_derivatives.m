clc;clear; close all;

% set seed
rng('default')

% define parameters
t0 = 0; tf = 1;

nt = 200; nsamples = 50;
fun_name = 'vanderpol';

% run simulation and get results
sim_details = run_simulation(t0,tf,nt,nsamples,fun_name);

split = [1,0];

% dfsm options
dfsm_options.ltype = 'LTI';
dfsm_options.ntype = 'RBF';
dfsm_options.nsamples = 500;
dfsm_options.sampling_type = 'KM';
dfsm_options.train_test_split = split;
dfsm_options.scale_flag = ~true;


dfsm = DFSM(sim_details,dfsm_options);

nstates = sim_details(1).nstates;
lm = dfsm.lin;

lm = lm';
x = lm(1:nstates,1:nstates);

eig(x)

%  
return



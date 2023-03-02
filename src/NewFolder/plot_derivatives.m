clc;clear; close all;

% set seed
rng('default')

% define parameters
t0 = 0; tf = 3;

nt = 200; nsamples = 10;
fun_name = 'two-link-robot';

% run simulation and get results
sim_details = run_simulation(t0,tf,nt,nsamples,fun_name);

split = [0.9,0.1];

% dfsm options
dfsm_options.ltype = 'LTI';
dfsm_options.ntype = 'NN';
dfsm_options.lsamples = nan;
dfsm_options.nsamples = 3000;
dfsm_options.sampling_type = 'KM';
dfsm_options.train_test_split = split;


dfsm = DFSM(sim_details,dfsm_options);


% 
% 
% % evaluate the difference between actual and polynomial approximate
% dx_act = vertcat(dXcell_act{:});
% dx_poly = vertcat(dXcell_poly{:});
% dx_diff = dx_act-dx_poly;
% 
% % plot 
% nplot = 5;
% hf = figure; hf.Color = 'w';
% hold on;
% plot(Tcell{nplot},dXcell_act{nplot},'k','linewidth',1)
% plot(Tcell{nplot},dXcell_poly{nplot},'r--','linewidth',1)
% 
% % histogram
% hf = figure;hf.Color = 'w';
% hold on;
% histogram(dx_diff(:,1))
% histogram(dx_diff(:,2))
%  
return



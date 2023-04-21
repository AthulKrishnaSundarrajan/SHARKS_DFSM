clc;clear; close all;

% set seed
rng(4357)

% define parameters
t0 = 0; tf = 3;

nt = 50; nsamples = 100;
fun_name = 'two-link-robot';

% run simulation and get results
sim_details = run_simulation(t0,tf,nt,nsamples,fun_name,1);

split = [1,0];
saveflag = false;

ind_test = randsample(1:nsamples,1);

time = sim_details(ind_test).time;
dx_poly = sim_details(ind_test).state_derivatives;
dx_act = sim_details(ind_test).dx_act;

nt = length(time);

ind = 1:nt;

ind_ = linspace(1,nt,100);

time_ = interp1(ind,time,ind_);
dx_poly_ = interp1(time,dx_poly,time_);
dx_act_ = interp1(time,dx_act,time_);


%---------------------------------------------------------
% hf = figure;
% hf.Color = 'w';
% hold on;
% hf.Position = [1314 469 570 413];
% commonFigureProperties
% 
% plot(time,dx_poly,'k','LineWidth',1.5)
% plot(time,dx_act,'r-','Linewidth',1.5)
% 
% xlabel('Time [s]');ylabel('$\dot{\xi}$')
% taskflag = 'axes'; commonFigureTasks;
% 
% legend('Polynomial approximation','','','','Actual');
% fontlegend = 22; nCol = 3;lcn = 'best';
% taskflag = 'legend';commonFigureTasks;
% 
% if saveflag
%     savename = 'control_inputs';
%     pathpdf = mfoldername(mfilename('fullpath'),fol_name);
%     filename = fullfile(pathpdf,savename);
%     str = strcat("export_fig '",filename,"' -pdf");
%     eval(str)
% 
% end
%--------------------------------------------------------
% K-means approach
[input_KM,dx_KM,output_KM,inputs,state_dx,outputs] = sample_data(sim_details,'KM',500);

% unique tol
[input_UTOL,dx_UTOL,~] = uniqueDataTol(inputs,state_dx,0.5);

% plot state space
hf = figure;
hf.Color = 'w';
hold on;

plot(inputs(:,3),inputs(:,4),'.','MarkerSize',8)
plot(input_KM(:,3),input_KM(:,4),'.','markersize',16)
%plot(input_UTOL(:,3),input_UTOL(:,4),'.','MarkerSize',16)

xlabel('\xi_1');ylabel('\xi_2')

legend({'All samples','Kmeans','Unique-Tol'},'NumColumns',3,'Location','northoutside')

return

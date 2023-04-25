clc;clear; close all;

% set seed
rng(4357)

% define parameters
t0 = 0; tf = 3;

nt = 90; nsamples = 100;
fun_name = 'two-link-robot';

% run simulation and get results
sim_details = run_simulation(t0,tf,nt,nsamples,fun_name,1);

split = [0.8,0.2];
saveflag = ~false;
fol_name = 'plot_sampling';

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
close all
% unique tol
%[input_UTOL,dx_UTOL,~] = uniqueDataTol(inputs,state_dx,0.5);

% plot state space
hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [1314 469 570 413];
commonFigureProperties

plot(inputs(:,3),inputs(:,4),'-','LineWidth',2)
plot(input_KM(:,3),input_KM(:,4),'.','markersize',20)
xlim([-2,2]);ylim([-2,2])
%plot(input_UTOL(:,3),input_UTOL(:,4),'.','MarkerSize',16)

xlabel('$\xi_1$');ylabel('$\xi_2$')
taskflag = 'axes';commonFigureTasks;

% export
if saveflag
    savename = 'states';
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end

%----------------------------------------------------------------------

hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [1314 469 570 413];
commonFigureProperties

plot(inputs(:,1),inputs(:,2),'-','Linewidth',2)
plot(input_KM(:,1),input_KM(:,2),'.','markersize',20)


xlabel('$u_1$');ylabel('$u_2$')
xlim([-1.5,1.5]);ylim([-2.5,2.5])
taskflag = 'axes';commonFigureTasks;

% export
if saveflag
    savename = 'control';
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end


%------------------------------------------------------
legend_entries = {'All samples','K-means'};

hf2 = copyfig(hf);
hf2.Position = [263.4000 417.8000 500.2000 70.4000];
h = findall(hf2,'type','line');
close(hf)

axis off; 
for i = 1:length(h)

    h(i).XData = NaN; %ignore warnings

end

% legend
legend(legend_entries);
fontlegend = 22; nCol = 3;lcn = 'best';
taskflag = 'legend';commonFigureTasks;


% export
if saveflag
    savename = 'legend_common';
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end

return

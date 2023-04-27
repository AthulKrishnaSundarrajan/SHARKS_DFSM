clc;clear; close all;

% set seed
rng(4357)

% define parameters
t0 = 0; tf = 3;

nt = 90; nsamples = 100;
fun_name = 'two-link-robot';

fac = 1; x0 = zeros(1,4);

% run simulation and get results
sim_details = run_simulation(t0,tf,nt,nsamples,fun_name,fac,x0);

split = [0.8,0.2];
saveflag = ~false;
fol_name = 'plot_sampling';

x0_array = zeros(0.8*nsamples,4);

for i = 1:0.8*nsamples

X = sim_details(i).states;

x0_array(i,:) = X(1,:);

end

C = materialColors;

cred = C.red(9,:);


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

plot(x0_array(:,3),x0_array(:,4),'k.','markersize',20)
plot(input_KM(:,5),input_KM(:,6),'.','color',cred,'markersize',20)
xlim([-1.5,1.5]);ylim([-1.5,1.5])
%plot(input_UTOL(:,3),input_UTOL(:,4),'.','MarkerSize',16)

xlabel('$\xi_3$');ylabel('$\xi_4$')
taskflag = 'axes';commonFigureTasks;

% export
if saveflag
    savename = 'states_x3x4';
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

plot(input_KM(:,1),input_KM(:,2),'r.','color',cred,'markersize',20)


xlabel('$u_1$');ylabel('$u_2$')
xlim([-1,1]);ylim([-1,1])
taskflag = 'axes';commonFigureTasks;

% export
if saveflag
    savename = 'control';
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end

%------------------------------------------------------------
% plot state space
hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [1314 469 570 413];
commonFigureProperties

plot(x0_array(:,1),x0_array(:,2),'k.','markersize',20)
plot(input_KM(:,3),input_KM(:,4),'.','color',cred,'markersize',20)
xlim([-1.2,1.2]);ylim([-1.2,1.2])
%plot(input_UTOL(:,3),input_UTOL(:,4),'.','MarkerSize',16)

xlabel('$\xi_1$');ylabel('$\xi_2$')
taskflag = 'axes';commonFigureTasks;

% export
if saveflag
    savename = 'states_x1x2';
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end


%----------------------------------------------------------------------


%------------------------------------------------------
legend_entries = {'Starting Points','K-means Samples'};

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

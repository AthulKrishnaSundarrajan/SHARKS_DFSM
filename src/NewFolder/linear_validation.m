clc; clear; close all;

% set seed
rng(4357)

% define parameters
t0 = 0; tf = 3;
fname = mfilename('fullpath');


fname = which(fname);

nt = 90; nsamples = 50;
fun_name = 'two-link-robot';
split = [0.7,0.3];

% run simulation and get results
nfac = 1;

sim_details_cell = cell(nfac,1);
dfsm_cell = cell(nfac,1);
inputs_cell = cell(nfac,1);


fac = 0.01; %logspace(-2,0,nfac);

dfsm_options.ltype = 'LTI';
dfsm_options.ntype = '';
dfsm_options.nsamples = nan;
dfsm_options.sampling_type = 'KM';
dfsm_options.train_test_split = split;
dfsm_options.scale_flag = ~true;


for i = 1:nfac
    sim_details_cell{i} = run_simulation(t0,tf,nt,nsamples,fun_name,fac(i));

    dfsm_cell{i} = DFSM(sim_details_cell{i},dfsm_options);

    [~,~,~,inputs,~,~] = sample_data(sim_details_cell{i},'KM',nan);

    inputs_cell{i} = inputs;

end


%[input_sampled,state_dx_sampled,output_sampled,inputs,state_dx,outputs] = sample_data(sim_details,'KM',nan);


saveflag = false;
fol_name = 'plots_linear_validation';
x_lim = [0,tf];


%-------------------------------------------------
% define symbolic variables
syms x1 x2 x3 x4 u1 u2 real

% derivative function
x1_d = ((sin(x3).*(9/4*cos(x3).*x1.^2+2*x2.^2) + 4/3*(u1-u2) - 3/2*cos(x3).*u2 )./ (31/36 + 9/4*sin(x3).^2) );
x2_d = (-( sin(x3).*(9/4*cos(x3).*x2.^2+7/2*x1.^2) - 7/3*u2 + 3/2*cos(x3).*(u1-u2) )./ (31/36 + 9/4*sin(x3).^2) );
x3_d = ( x2-x1 ); 
x4_d = ( x1 );

% outputs
dx = [x1_d;x2_d;x3_d;x4_d];

% inputs
x = [x1;x2;x3;x4];
u = [u1;u2];

% inputs
I = [u;x];

% evaluate jacobian
L_taylor = jacobian(dx,I);

% substitute value
L_taylor = double(subs(L_taylor,I,zeros(6,1)));

%-------------------------------------------------

% create a dfsm model with the taylor series expansion as the linear model
dfsm_taylor = dfsm_cell{end};
dfsm_taylor.deriv.AB = L_taylor';


mse_dx_taylor = zeros(nfac,nsamples,length(x));
mse_dx_linfit = zeros(nfac,nsamples,length(x));
lin_model = cell(nfac,1);

for i = 1:nfac

    sim_details = sim_details_cell{i};
    dfsm_ = dfsm_cell{i};
    lin_model{i} = dfsm_.deriv.AB;

    [~,~,dx_linfit] = test_dfsm(dfsm_cell{i},sim_details,1:nsamples,false,false);
    [~,~,dx_taylor] = test_dfsm(dfsm_taylor,sim_details,1:nsamples,false,false);

    for n = 1:nsamples

        mse_dx_taylor(i,n,:) = calculate_mse(dx_taylor{n,1},dx_taylor{n,2});
        mse_dx_linfit(i,n,:) = calculate_mse(dx_linfit{n,1},dx_linfit{n,2});


    end

end

mat_name = 'mse-data.mat';
save(mat_name,'fac','nsamples',"mse_dx_taylor",'mse_dx_linfit',"lin_model","L_taylor","inputs_cell");



%--------------------------------------------------

fac_ind = nfac;
dfsm = dfsm_cell{fac_ind};
sim_details = sim_details_cell{nfac};

% extract linear model
L_dfsm = (dfsm.deriv.AB)';
eig3 = eig(L_dfsm(:,3:end));



% generate a random index to test
ind_test = randsample(1:nsamples,1);

% extract values for this index
time = sim_details(ind_test).time;
states = sim_details(ind_test).states;
controls = sim_details(ind_test).controls;
state_derivatives = sim_details(ind_test).state_derivatives;

% test dfsm for the index
[dfsm,Xcell_dfsm,dx_dfsm] = test_dfsm(dfsm,sim_details(ind_test),1,false,true);
[dfsm_taylor,Xcell_taylor,dx_taylor] = test_dfsm(dfsm_taylor,sim_details(ind_test),1,false,true);

% extract state values
X_dfsm = Xcell_dfsm{1,2};
X_taylor= Xcell_taylor{1,2};
%--------------------------------------------
close all;
hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [1314 469 570 413];
commonFigureProperties

plot(time,controls(:,1),'LineWidth',1.5)
plot(time,controls(:,2),'Linewidth',1.5)

xlabel('Time [s]');ylabel('$u$')
taskflag = 'axes'; commonFigureTasks;

legend('$u_1$','$u_2$');
fontlegend = 22; nCol = 3;lcn = 'best';
taskflag = 'legend';commonFigureTasks;

if saveflag
    savename = 'control_inputs';
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)

end


%--------------------------------------------

% 
legend_entries = {'$\hat{f}_{\textrm{linfit}}$','$\hat{f}_{\textrm{taylor}}$','$f_{\textrm{actual}}$'};

%-------------------------------

for idx = 1:4
    hf = figure;
    hf.Color = 'w';
    hold on;
    hf.Position = [1314 469 570 413];
    commonFigureProperties
    
    y_label = ['$\xi_',num2str(idx),'$'];
    
    
    plot(time,X_dfsm(:,idx),'linewidth',1.5)
    plot(time,X_taylor(:,idx),'linewidth',1.5)
    plot(time,states(:,idx),'linewidth',1.5)

    xlim(x_lim)

    taskflag = 'axes'; commonFigureTasks;
    
    xlabel('Time [s]');ylabel(y_label)

    if saveflag
        savename = ['state_traj_',num2str(idx)];
        pathpdf = mfoldername(mfilename('fullpath'),fol_name);
        filename = fullfile(pathpdf,savename);
        str = strcat("export_fig '",filename,"' -pdf");
        eval(str)

    end
end
%------------------------------------------

dx_dfsm = dx_dfsm{1,2};
dx_taylor = dx_taylor{1,2};

for idx = 1:4

    hf = figure;
    hf.Color = 'w';
    hold on;
    hf.Position = [1314 469 570 413];
    commonFigureProperties
    
    y_label = ['$\dot{\xi}_',num2str(idx),'$'];

    plot(time,dx_dfsm(:,idx),'linewidth',1.5)
    plot(time,dx_taylor(:,idx),'linewidth',1.5)
    plot(time,state_derivatives(:,idx),'linewidth',1.5)
    xlim(x_lim)

    xlabel('Time [s]');ylabel(y_label)
    
    taskflag = 'axes'; commonFigureTasks;

    if saveflag
        savename = ['dx_traj_',num2str(idx)];
        pathpdf = mfoldername(mfilename('fullpath'),fol_name);
        filename = fullfile(pathpdf,savename);
        str = strcat("export_fig '",filename,"' -pdf");
        eval(str)

    end

   

end

%----------------------
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
    savename = 'legend_common_traj';
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end


%---------------------------------------------
%%
%[input_sampled,state_dx_sampled,output_sampled,inputs,state_dx,outputs] = sample_data(sim_details,'KM',100);



mse_dx_taylor = zeros(nfac,nsamples);
mse_dx_linfit = zeros(nfac,nsamples);



dx_dfsm = vertcat(dx_dfsm{:,2});
dx_taylor = vertcat(dx_taylor{:,2});
dx_actual = vertcat(sim_details(:).state_derivatives);



diff_linfit = dx_dfsm - dx_actual;
diff_taylor = dx_taylor - dx_actual;

nbins = 70;
x_lim = {[-0.1,0.1],[-0.1,0.1],[-1e-7,1e-7],[-1e-7,1e-7]};

for idx = 1:2

hf = figure;
hf.Color = 'w';
hold on;
hf.Position = [1314 469 570 413];
commonFigureProperties

x_label = ['$\dot{\xi}_',num2str(idx),'$'];

histogram(diff_linfit(:,idx),'FaceAlpha',0.7,'NumBins',nbins,'Normalization','count')
histogram(diff_taylor(:,idx),'FaceAlpha',0.7,'NumBins',nbins,'Normalization','count')

xlim(x_lim{idx}); xlabel(x_label);
ylabel('Count')

taskflag = 'axes'; commonFigureTasks;

if saveflag
        savename = ['dx_hist_',num2str(idx)];
        pathpdf = mfoldername(mfilename('fullpath'),fol_name);
        filename = fullfile(pathpdf,savename);
        str = strcat("export_fig '",filename,"' -pdf");
        eval(str)

end



end
%----------------------
hf2 = copyfig(hf);
hf2.Position = [263.4000 417.8000 500.2000 70.4000];
h = findall(hf2,'type','histogram');
close(hf)
axis off; 
for i = 1:length(h)

    h(i).Data = NaN; %ignore warnings

end

% legend
legend({'$e_{\textrm{linfit}}$','$e_{\textrm{taylor}}$'});
fontlegend = 22; nCol = 2;lcn = 'best';
taskflag = 'legend';commonFigureTasks;


% export
if saveflag
    savename = 'legend_common_hist';
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end

return
clc; clear; close all;

% set seed
rng(4357)

% define parameters
t0 = 0; t_f = 3;
fname = mfilename('fullpath');


fname = which(fname);

nt = 90; nsamples = 100;
fun_name = 'two-link-robot';
split = [1,0];

% run simulation and get results

fac = 0.1; %logspace(-2,0,nfac);
x0 = [0,0,0.5,0];

dfsm_options.ltype = 'LTI';
dfsm_options.ntype = '';
dfsm_options.nsamples = nan;
dfsm_options.sampling_type = 'KM';
dfsm_options.train_test_split = split;
dfsm_options.scale_flag = ~true;


sim_details = run_simulation(t0,t_f,nt,nsamples,fun_name,fac,x0);

train_simulations_ind = 1:floor(nsamples*0.8);
train_simulations = sim_details(train_simulations_ind);

test_simulations_ind = floor(nsamples*0.8+1):nsamples;
test_simulations = sim_details(test_simulations_ind);

dfsm = DFSM(train_simulations,dfsm_options);

[~,~,~,inputs,~,~] = sample_data(sim_details,'KM',nan);

saveflag = false;
fol_name = 'plots_linear_validation';
x_lim = [0,t_f];


%-------------------------------------------------
% define symbolic variables
syms x1 x2 x3 x4 u1 u2 X1 X2 X3 X4 U1 U2 real

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

x0 = [0;0;0.5;0];
u0 = [0;0];

I0 = [u0;x0];


% inputs
I = [u;x];
I_ = [U1;U2;X1;X2;X3;X4];

% evaluate jacobian
L_taylor = jacobian(dx,I);

% substitute value
L_taylor = double(subs(L_taylor,I,I0));

%--------------------------------------------------------------------------

% create a dfsm model with the taylor series expansion as the linear model

mse_dx_taylor = zeros(0.2*nsamples,length(x));
mse_dx_linfit = zeros(0.2*nsamples,length(x));

A_taylor = L_taylor(:,3:end);
B_taylor = L_taylor(:,1:2);


L_linfit = (dfsm.deriv.AB)';

[~,~,dx_linfit] = test_dfsm(dfsm,test_simulations,1:0.2*nsamples,false,false);
%[~,~,dx_taylor] = test_dfsm(dfsm_taylor,test_simulations,1:0.3*nsamples,false,false);

dx_taylor = test_taylor(A_taylor,B_taylor,test_simulations,x0,u0);

for n = 1:0.2*nsamples

    mse_dx_taylor(n,:) = calculate_mse(dx_taylor{n,1},dx_taylor{n,2});
    mse_dx_linfit(n,:) = calculate_mse(dx_linfit{n,1},dx_linfit{n,2});


end

A_linfit = L_linfit(:,3:end);
B_linfit = L_linfit(:,1:2);

ind_test = randsample(1:0.2*nsamples,1);

dx_act = dx_taylor{ind_test,1};
dx_taylor1 = dx_taylor{ind_test,2};
dx_linfit1 = dx_linfit{ind_test,2};


time = test_simulations(ind_test).time;
%--------------------------------------------------------------------------
close all;
for idx = 1:4
    hf = figure;
    hf.Color = 'w'; hold on;
    commonFigureProperties
    
    
    plot(time,dx_linfit1(:,idx),'linewidth',1.5)
    plot(time,dx_taylor1(:,idx),'linewidth',1.5)
    plot(time,dx_act(:,idx),'linewidth',1.5)

    y_label = ['$\dot{\xi}_',num2str(idx),'$'];

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

legend_entries = {'$\hat{f}_{\textrm{linfit}}$','$\hat{f}_{\textrm{taylor}}$','$f_{\textrm{actual}}$'};

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

%-------------------------------------------------------------------------

C_ = [0 0 1 0;
    0 0 0 1];
D = zeros(2);

% A_linfit(abs(A_linfit)<1e-2) = 0;
% B_linfit(abs(B_linfit)<1e-2) = 0;
sys_linfit = ss(A_linfit,B_linfit,C_,D);
sys_taylor = ss(A_taylor,B_taylor,C_,D);

tf_linfit = tf(sys_linfit);tf_taylor = tf(sys_taylor);


W = logspace(-1,3,1e4);

[MAG1,PHASE1] = bode(sys_linfit,W);

% extract 
x1u1_linfit = squeeze(MAG1(1,1,:));
x1u2_linfit = squeeze(MAG1(1,2,:));
x2u1_linfit = squeeze(MAG1(2,1,:));
x2u2_linfit = squeeze(MAG1(2,2,:));

% evaluate frequency domain result
[MAG2,PHASE2] = bode(sys_taylor,W);

x1u1_taylor = squeeze(MAG2(1,1,:));
x1u2_taylor = squeeze(MAG2(1,2,:));
x2u1_taylor = squeeze(MAG2(2,1,:));
x2u2_taylor = squeeze(MAG2(2,2,:));

%--------------------------------------------------------------------------
% x1u1
%--------------------------------------------------------------------------
hf = figure;
hf.Color = 'w'; hold on;
commonFigureProperties
C = materialColors;


plot(W,mag2db(x1u1_linfit),'linewidth',2)
plot(W,mag2db(x1u1_taylor),'linewidth',2)

% x and y label
xlabel('$\omega$ [rad/s]');xlim([W(1),W(end)])
ylabel('Magnitude [dB]'); 

% axes
taskflag = 'axes'; commonFigureTasks;
ha.XScale = 'log';ha.YScale = 'linear';

if saveflag
        savename = ['freq_x1u1_',num2str(idx)];
        pathpdf = mfoldername(mfilename('fullpath'),fol_name);
        filename = fullfile(pathpdf,savename);
        str = strcat("export_fig '",filename,"' -pdf");
        eval(str)

end

%--------------------------------------------------------------------------
% x1u2
%--------------------------------------------------------------------------
hf = figure;
hf.Color = 'w'; hold on;
commonFigureProperties
C = materialColors;


plot(W,mag2db(x1u2_linfit),'linewidth',2)
plot(W,mag2db(x1u2_taylor),'linewidth',2)

% x and y label
xlabel('$\omega$ [rad/s]');xlim([W(1),W(end)])
ylabel('Magnitude [dB]'); 

% axes
taskflag = 'axes'; commonFigureTasks;
ha.XScale = 'log';ha.YScale = 'linear';

if saveflag
        savename = ['freq_x1u2_',num2str(idx)];
        pathpdf = mfoldername(mfilename('fullpath'),fol_name);
        filename = fullfile(pathpdf,savename);
        str = strcat("export_fig '",filename,"' -pdf");
        eval(str)

end

%--------------------------------------------------------------------------
% x2u1
%--------------------------------------------------------------------------
hf = figure;
hf.Color = 'w'; hold on;
commonFigureProperties
C = materialColors;


plot(W,mag2db(x2u1_linfit),'linewidth',2)
plot(W,mag2db(x2u1_taylor),'linewidth',2)

% x and y label
xlabel('$\omega$ [rad/s]');xlim([W(1),W(end)])
ylabel('Magnitude [dB]'); 

% axes
taskflag = 'axes'; commonFigureTasks;
ha.XScale = 'log';ha.YScale = 'linear';


if saveflag
        savename = ['freq_x2u1_',num2str(idx)];
        pathpdf = mfoldername(mfilename('fullpath'),fol_name);
        filename = fullfile(pathpdf,savename);
        str = strcat("export_fig '",filename,"' -pdf");
        eval(str)

end


%--------------------------------------------------------------------------
% x2u2
%--------------------------------------------------------------------------
hf = figure;
hf.Color = 'w'; hold on;
commonFigureProperties
C = materialColors;


plot(W,mag2db(x2u2_linfit),'linewidth',2)
plot(W,mag2db(x2u2_taylor),'linewidth',2)

% x and y label
xlabel('$\omega$ [rad/s]');xlim([W(1),W(end)])
ylabel('Magnitude [dB]'); 

% axes
taskflag = 'axes'; commonFigureTasks;
ha.XScale = 'log';ha.YScale = 'linear';

if saveflag
        savename = ['freq_x2u2_',num2str(idx)];
        pathpdf = mfoldername(mfilename('fullpath'),fol_name);
        filename = fullfile(pathpdf,savename);
        str = strcat("export_fig '",filename,"' -pdf");
        eval(str)

end



%--------------------------------------------------------------------------

legend_entries = {'$\hat{f}_{\textrm{linfit}}$','$\hat{f}_{\textrm{taylor}}$'};
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
    savename = 'legend_common_freqdom';
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end
%-------------------------------------------------------------------------------

% mat_name = 'mse-data.mat';
% save(mat_name,'fac','nsamples',"mse_dx_taylor",'mse_dx_linfit',"lin_model","L_taylor","inputs_cell");
function dx_taylor = test_taylor(A,B,test_simulations,x0,u0)

    ntest = length(test_simulations);

    dx_taylor = cell(ntest,2);

    for i = 1:ntest

        states = test_simulations(i).states;
        controls = test_simulations(i).controls;
        state_derivatives = test_simulations(i).state_derivatives;

        dx = (states-x0' )*A' + (controls-u0')*B';

        dx_taylor{i,1} = state_derivatives;
        dx_taylor{i,2} = dx;
        

    end
    
end
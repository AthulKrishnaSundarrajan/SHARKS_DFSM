clc; clear; close all;

% load the data
load('TLR_state_samples_LHC.mat');

% plot state space samples
hf = figure;
hf.Color = 'w';
hold on;
commonFigureProperties;

C = materialColors;
c_all = C.blue(10,:);
c_samp = C.red(10,:);

markersize = 18;
plot(X_sampled(:,1),X_sampled(:,2),'.','markersize',2,'color',c_all);
plot(input_sampled(:,3),input_sampled(:,4),'.','markersize',markersize,'color',c_samp)
xlim([-1.5,1.5]);ylim([-1.5,1.5]);
xlabel('$\xi_1$'); ylabel('$\xi_2$')

%legend('')

taskflag = 'axes';commonFigureTasks

% plot control inputs
hf = figure;
hf.Color = 'w';
hold on;
commonFigureProperties

plot(input_sampled(:,1),input_sampled(:,2),'.','markersize',markersize,'color',c_samp)
xlim([-1,1]);ylim([-1,1]);
xlabel('$u_1$');ylabel('$u_2$')
taskflag = 'axes';commonFigureTasks
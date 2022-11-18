clc; clear; close all;

% load results
results = load('lstm_weis_results.mat');

% extract results
time = results.time;time(end) = [];nt = length(time);
time_test = results.time_test;
time_train = results.time_train;

indTrain = 1: floor(0.4*nt);ntrain = length(indTrain);
indTest = indTrain(end)+1:nt; ntest = length(indTest);

% states
states = results.states;states(end,:) = [];
controls = results.controls;controls(end,:) = [];
state_derivatives  = results.state_derivatives;state_derivatives(end,:) = [];

% training and testing data
XTest = results.XTest;
XTrain = results.XTrain;

TTest = results.TTest;
TTrain = results.TTrain;

dx_lstm = results.dx_lstm;
muT = results.muT; muX = results.muX;
sigmaT = results.sigmaT; sigmaX = results.sigmaX;

dx_lstm = dx_lstm.*sigmaT;
dx_lstm = dx_lstm + muT;

%% inputs
% plot
hf = figure;
hf.Color = 'w';

subplot(4,2,1)
hold on;
plot(time_train,controls(indTrain,1),'k')
plot(time_test,controls(indTest,1),'r')
xlim([time(1),time(end)])
%legend('Training','Validation')
title('Wind Speed [m/s]')

subplot(4,2,3)
hold on;
plot(time_train,controls(indTrain,2)/1000,'k')
plot(time_test,controls(indTest,2)/1000,'r')
xlim([time(1),time(end)])
title('Gen Torque [MWm]')

subplot(4,2,5)
hold on;
plot(time_train,controls(indTrain,3),'k')
plot(time_test,controls(indTest,3),'r')
xlim([time(1),time(end)])
title('Blade Pitch [deg]')

subplot(4,2,7)
hold on;
plot(time_train,states(indTrain,1),'k')
plot(time_test,states(indTest,1),'r')
xlim([time(1),time(end)])
title('Platform Pitch [deg]')


subplot(4,2,2)
hold on;
plot(time_train,states(indTrain,2),'k')
plot(time_test,states(indTest,2),'r')
%legend('Training','Validation','NumColumns',2)
xlim([time(1),time(end)])
title('Gen Speed [rpm]')

subplot(4,2,4)
hold on;
plot(time_train,states(indTrain,3),'k')
plot(time_test,states(indTest,3),'r')
%legend('Training','Validation','NumColumns',2)
xlim([time(1),time(end)])
title('dPtfmPitch')

subplot(4,2,6)
hold on;
plot(time_train,states(indTrain,4),'k')
plot(time_test,states(indTest,4),'r')
legend('Training','Validation','NumColumns',2)
xlim([time(1),time(end)])
title('dGenSpeed')

sgtitle('Inputs')
%% outputs

hf = figure;
hf.Color = 'w';

ind = 1;

subplot(2,2,ind)
hold on;
plot(time_train,state_derivatives(indTrain,ind),'k')
plot(time_test,state_derivatives(indTest,ind),'r')
xlabel('Time'); xlim([time(1),time(end)])
title('dPtfmPitch')

ind = ind + 1;


subplot(2,2,ind)
hold on;
plot(time_train,state_derivatives(indTrain,ind),'k')
plot(time_test,state_derivatives(indTest,ind),'r')
xlabel('Time'); xlim([time(1),time(end)])
title('dGenSpeed')

ind = ind + 1;


subplot(2,2,ind)
hold on;
plot(time_train,state_derivatives(indTrain,ind),'k')
plot(time_test,state_derivatives(indTest,ind),'r')
xlabel('Time'); xlim([time(1),time(end)])
title('ddPtfmPitch')


ind = ind + 1;


subplot(2,2,ind)
hold on;
plot(time_train,state_derivatives(indTrain,ind),'k')
plot(time_test,state_derivatives(indTest,ind),'r')
xlabel('Time'); xlim([time(1),time(end)])
title('ddGenSpeed')

sgtitle('Outputs')

%% results

hf = figure;
hf.Color = 'w';
ind = 1;

subplot(2,2,ind)
hold on;
%plot(time_test,TTest(ind,:),'k')
plot(time_test,state_derivatives(indTest,ind),'k')
plot(time_test,dx_lstm(ind,:),'r--')
title('dPtfmPitch')
xlim([time_test(1),time_test(end)])
%legend('original','DFSM LSTM')


ind = ind + 1;

subplot(2,2,ind)
hold on;
%plot(time_test,TTest(ind,:),'k')
plot(time_test,state_derivatives(indTest,ind),'k')
plot(time_test,dx_lstm(ind,:),'r--')
title('dGenSpeed')
xlim([time_test(1),time_test(end)])

ind = ind + 1;

subplot(2,2,ind)
hold on;
%plot(time_test,TTest(ind,:),'k')
plot(time_test,state_derivatives(indTest,ind),'k')
plot(time_test,dx_lstm(ind,:),'r--')
title('ddPtfmPitch')
xlim([time_test(1),time_test(end)])


ind = ind + 1;

subplot(2,2,ind)
hold on;
%plot(time_test,TTest(ind,:),'k')
plot(time_test,state_derivatives(indTest,ind),'k')
plot(time_test,dx_lstm(ind,:),'r--')
title('ddGenSpeed')
xlim([time_test(1),time_test(end)])
legend('Original','LSTM','NumColumns',2)

sgtitle('Results')
%% 

clc; clear; close all;

% load results
results = load('lstm_weis_results_60red.mat');

% extract results
time = results.time;%time(end) = [];
nt = length(time);
time_test = results.time_test;nsplit = results.nsplit;
time_train = results.time_train;

indTrain = 1: floor(nsplit/100*nt);ntrain = length(indTrain);
indTest = indTrain(end)+1:nt-1; ntest = length(indTest);

% states
states = results.states;%states(end,:) = [];
controls = results.controls;%controls(end,:) = [];
state_derivatives  = results.state_derivatives;%state_derivatives(end,:) = [];

% training and testing data
XTest = results.XTest;
XTrain = results.XTrain;

TTest = results.TTest;
TTrain = results.TTrain;

dx_lstm = results.dx_lstm;
muT = results.muT; muX = results.muX;
sigmaT = results.sigmaT; sigmaX = results.sigmaX;

T_OF = results.T_OF; X_OF = results.X_OF;
T_dfsm = results.T_dfsm;X_dfsm = results.X_dfsm;



%% inputs
% plot
hf = figure;
hf.Color = 'w';

subplot(4,2,1)
hold on;
plot(time_train,controls(indTrain,1),'k')
plot(time_test,controls(indTest,1),'r')
xlim([time(1),time(end)])
xlabel('Time')
%legend('Training','Validation')
title('Wind Speed [m/s]')

subplot(4,2,3)
hold on;
plot(time_train,controls(indTrain,2)/1000,'k')
plot(time_test,controls(indTest,2)/1000,'r')
xlim([time(1),time(end)])
xlabel('Time')
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
xlabel('Time')
title('Platform Pitch [deg]')


subplot(4,2,2)
hold on;
plot(time_train,states(indTrain,2),'k')
plot(time_test,states(indTest,2),'r')
%legend('Training','Validation','NumColumns',2)
xlim([time(1),time(end)])
xlabel('Time')
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
%legend('Training','Validation','NumColumns',2)
xlim([time(1),time(end)])
xlabel('Time')
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

%% common legend

hf2 = figure;
hf2.Color = 'w';
hold on;

%close(hf);

plot(1:10,1:10,'r')
plot(1:5,1:5,'k')
h = findall(hf2,'type','line');


legend([h(1) h(2)],{'Training','Validation'},'fontsize',18,'NumColumns',2,'Location','northoutside');
fonttick = 25; nCol = 3;lcn = 'northoutside';
%taskflag = 'legend';commonFigureTasks;

axis off; 
for i = 1:length(h)

    h(i).XData = NaN; %ignore warnings

end

hf2.Position = [1387 834 478 57];

%% results

hf = figure;
hf.Color = 'w';
ind = 1;

subplot(2,2,ind)
hold on;
plot(T_OF,X_OF(:,ind),'k')
plot(T_dfsm,X_dfsm(:,ind),'r-')
title('PtfmPitch'); xlabel('Time [s]')
xlim([time_test(1),time_test(end)])
%legend('original','DFSM LSTM')


ind = ind + 1;

subplot(2,2,ind)
hold on;
%plot(time_test,TTest(ind,:),'k')
plot(T_OF,X_OF(:,ind),'k')
plot(T_dfsm,X_dfsm(:,ind),'r-')
title('GenSpeed'); xlabel('Time [s]')
xlim([time_test(1),time_test(end)])

ind = ind + 1;

subplot(2,2,ind)
hold on;
%plot(time_test,TTest(ind,:),'k')
plot(T_OF,X_OF(:,ind),'k')
plot(T_dfsm,X_dfsm(:,ind),'r-')
title('dPtfmPitch'); xlabel('Time [s]')
xlim([time_test(1),time_test(end)])


ind = ind + 1;

subplot(2,2,ind)
hold on;
%plot(time_test,TTest(ind,:),'k')
plot(T_OF,X_OF(:,ind),'k')
plot(T_dfsm,X_dfsm(:,ind),'r-')
title('dGenSpeed'); xlabel('Time [s]')
xlim([time_test(1),time_test(end)])
%legend('Original','LSTM','NumColumns',2)

sgtitle('State Trajectories')
%% dx 
hf = figure;
hf.Color = 'w';
sgtitle('Predicted State Derivatives')

ind = 1;

subplot(2,2,ind)
hold on;
plot(time_test,state_derivatives(indTest,ind),'k')
plot(time_test,dx_lstm(ind,:),'r-')
title('dPtfmPitch'); xlabel('Time [s]')
%legend('original','DFSM LSTM')


ind = ind + 1;

subplot(2,2,ind)
hold on;
plot(time_test,state_derivatives(indTest,ind),'k')
plot(time_test,dx_lstm(ind,:),'r-')
title('dGenSpeed'); xlabel('Time [s]')
%legend('original','DFSM LSTM')


ind = ind + 1;

subplot(2,2,ind)
hold on;
plot(time_test,state_derivatives(indTest,ind),'k')
plot(time_test,dx_lstm(ind,:),'r-')
title('ddPtfmPitch'); xlabel('Time [s]')
%legend('original','DFSM LSTM')


ind = ind + 1;

subplot(2,2,ind)
hold on;
plot(time_test,state_derivatives(indTest,ind),'k')
plot(time_test,dx_lstm(ind,:),'r-')
title('ddGenSpeed'); xlabel('Time [s]')
%legend('original','DFSM LSTM')

%%
hf2 = figure;
hf2.Color = 'w';
hold on;

%close(hf);

plot(1:10,1:10,'r')
plot(1:5,1:5,'k')
h = findall(hf2,'type','line');


legend([h(1) h(2)],{'OpenFAST','DFSM'},'fontsize',18,'NumColumns',2,'Location','northoutside');
fonttick = 25; nCol = 3;lcn = 'northoutside';
%taskflag = 'legend';commonFigureTasks;

axis off; 
for i = 1:length(h)

    h(i).XData = NaN; %ignore warnings

end

hf2.Position = [1387 834 478 57];
%% 

err_curve = X_OF - X_dfsm;

hf = figure;
hf.Color = 'w';
ind = 1;

subplot(2,2,ind)
hold on;
histogram(err_curve(:,ind))
title('PtfmPitch')
%xlim([time_test(1),time_test(end)])


ind = ind + 1;

subplot(2,2,ind)
hold on;
histogram(err_curve(:,ind))
title('GenSpeed')
%xlim([time_test(1),time_test(end)])

ind = ind + 1;

subplot(2,2,ind)
hold on;
histogram(err_curve(:,ind))
title('dPtfmPitch')
%xlim([time_test(1),time_test(end)])


ind = ind + 1;

subplot(2,2,ind)
hold on;
histogram(err_curve(:,ind))
title('dGenSpeed')

sgtitle('$x_{\textrm{dfsm}}-x_{\textrm{OF}}$','interpreter','latex','FontSize',15)

%%
err_curve = dx_lstm' - state_derivatives(indTest,:);

hf = figure;
hf.Color = 'w';
ind = 1;

subplot(2,2,ind)
hold on;
histogram(err_curve(:,ind))
title('dPtfmPitch')
%xlim([time_test(1),time_test(end)])


ind = ind + 1;

subplot(2,2,ind)
hold on;
histogram(err_curve(:,ind))
title('dGenSpeed')
%xlim([time_test(1),time_test(end)])

ind = ind + 1;

subplot(2,2,ind)
hold on;
histogram(err_curve(:,ind))
title('ddPtfmPitch')
%xlim([time_test(1),time_test(end)])


ind = ind + 1;

subplot(2,2,ind)
hold on;
histogram(err_curve(:,ind))
title('ddGenSpeed')

sgtitle('$dx_{\textrm{dfsm}}-dx_{\textrm{OF}}$','interpreter','latex','FontSize',15)

%% plot power
hf = figure;
hf.Color = 'w';
hold on;

controls_fun = @(t) interp1(time,controls,t,"pchip");

controls_ = controls_fun(T_OF);

p_dfsm = X_dfsm(:,2).*controls_(:,2)/1e4;
p_OF = X_OF(:,2).*controls_(:,2)/1e4;
p_diff = p_OF - p_dfsm;

plot(T_dfsm,p_dfsm,'r','linewidth',1)
plot(T_OF,p_OF,'k','linewidth',1)

xlim([T_OF(1),T_OF(end)])
xlabel('Time','fontsize',15); ylabel('Generator Power [MW]','fontsize',15)
legend('DFSM','OpenFAST')

return

function dx = dx_function(t,x,dx_fun)

    dx = dx_fun(t);
    dx = dx';

end

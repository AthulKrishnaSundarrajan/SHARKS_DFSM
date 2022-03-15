clc; clear; close all;

data = loadData_IEA_w_TMD;

data = approximateStateDerivatives(data);

%
% hf = figure; hold on
% hf.Color = 'w';
% for iCase = 1:1 % nLinCases
%
%     plot(data(iCase).time,data(iCase).states(:,1),'.','linewidth',1,...
%         'DisplayName',strcat("DLC",string(iCase)))
%
%     t_more = linspace(min(t),max(t),1e6);
%     x = ppval(pp_x_{iCase},t_more);
%     plot(t_more,x(1,:),'linewidth',1)
%
% end
% legend()
% xlim([200 210])
% ha = gca;
% ha.FontSize = 14;
% ha.LineWidth = 1;
% xlabel('time [s]')
% ylabel('')
%
% hf = figure; hold on
% hf.Color = 'w';
% for iCase = 1:1 % nLinCases
%
%     t = data(iCase).time;
%
%     plot(data(iCase).time,data(iCase).state_derivatives(:,2),'linewidth',1,...
%         'DisplayName',strcat("DLC",string(iCase)))
%
%     t_more = linspace(min(t),max(t),1e6);
%     Dx = ppval(pp_Dx_{iCase},t_more);
%     plot(t_more,Dx(2,:),'linewidth',1)
%
% end
% legend()
% xlim([200 210])
% ha = gca;
% ha.FontSize = 14;
% ha.LineWidth = 1;
% xlabel('time [s]')
% ylabel('')
%
%
%
%
%
% hf = figure; hold on
% hf.Color = 'w';
% for iCase = 1:1 % nLinCases
%
%     plot(data(iCase).time,data(iCase).states(:,2),'.','linewidth',1,...
%         'DisplayName',strcat("DLC",string(iCase)))
%
%     t_more = linspace(min(t),max(t),1e6);
%     x = ppval(pp_x_{iCase},t_more);
%     plot(t_more,x(2,:),'linewidth',1)
%
% end
% legend()
% xlim([200 210])
% ha = gca;
% ha.FontSize = 14;
% ha.LineWidth = 1;
% xlabel('time [s]')
% ylabel('')
%
% hf = figure; hold on
% hf.Color = 'w';
% for iCase = 1:1 % nLinCases
%
%     t = data(iCase).time;
%
%     plot(data(iCase).time,data(iCase).state_derivatives(:,1),'linewidth',1,...
%         'DisplayName',strcat("DLC",string(iCase)))
%
%     t_more = linspace(min(t),max(t),1e6);
%     Dx = ppval(pp_Dx_{iCase},t_more);
%     plot(t_more,Dx(2,:),'linewidth',1)
%
% end
% legend()
% xlim([200 210])
% ha = gca;
% ha.FontSize = 14;
% ha.LineWidth = 1;
% xlabel('time [s]')
% ylabel('')

% combine DLCs
Time = vertcat(data(:).time);
States = vertcat(data(:).states);
Inputs = vertcat(data(:).inputs);
State_derivatives = vertcat(data(:).state_derivatives);

% data for fitting
X = [States,Inputs];
Y = State_derivatives;


% find all wind speeds in range
p = 0.001; m = 4;
w = Inputs(:,1);
I = (m+p >= w) & (m-p <= w);
XI = X(I,:);
YI = Y(I,:);



Y1 = YI(:,1);
Y2 = YI(:,2);
Y1X = [XI,Y1];

% X = [ones(size(x1)) x1 x2 x1.*x2];

XI_ = [ones(size(Y1)) XI];
b1 = regress(Y1,XI_); % Removes NaN data
b2 = regress(Y2,XI_); % Removes NaN data

Y1_lin = XI_*b1;
Y2_lin = XI_*b2;

figure; hold on
plot(Y1_lin,Y1,'.')
plot(Y2_lin,Y2,'.')

return

z1 = iddata(Y,X,Time(2)-Time(1));
[sys,x0] = ssest(z1,1,'DisturbanceModel','none');

Y_ = lsim(sys,X,Time,x0);

figure; hold on
% plot(Y_,States,'.')
plot(Time,States,'.')
plot(Time,Y_,'.')

%

%
% X_ = tonndata(X,false,false);
% Y_ = tonndata(Y,false,false);
% net = timedelaynet(1:2,10);
% [Xs,Xi,Ai,Ts] = preparets(net,X_,Y_);
% net = train(net,Xs,Ts,Xi,Ai);
%
% [Y_,Xf,Af] = net(Xs,Xi,Ai);
% perf = perform(net,Ts,Y_);
%
% [netc,Xic,Aic] = closeloop(net,Xf,Af);
% view(netc)
%
% y2 = netc(X_,Xic,Aic);
%
% for k = 1:length(y2)
%     Y2(k,:) = y2{k};
% end
%
% figure; hold on
% plot(Y,Y2,'.')
%
% return
%
% [net,tr,xi,ai,t] = nndynamic(X,Y)



%
% Y_fit = net(X')';
% plot(Y,Y_fit,'.')
%
% netc = closeloop(net);
%
% Yc = netc(X');
% plot(Y,Yc')

return



Y1 = Y(:,1);
Y2 = Y(:,2);

[X_,Y_,maxX,maxY_] = scaleData(X,Y);

return


ny = size(Y,2);

for k = 1:ny


    %
    [X_,Y_,maxX,maxY_] = scaleData(X,Y(:,k));
    [X_,Y_] = uniqueDataTol(X_,Y_,0.05);
%     X_ = unscaleData(X_,maxX);
%     Y_ = unscaleData(Y_,maxY);

%     net = nntrain(X_,Y_);
%     net = newrb(X_',Y_',0,1,400);
%     net = newrbe(X_',Y_',100);
%     net = newgrnn(X_',Y_',0.1);
    [trainedModel, validationRMSE] = trainRegressionModel(X_,Y_);
    net = trainedModel.predictFcn2;
    models_{k} = net;

    figure
    try
        Y_fit = net(X_);
        plot(Y_,Y_fit,'.')
    catch
        Y_fit = net(X_');
        plot(Y_,Y_fit,'.')
    end


    [X_,Y_,maxX,maxY_] = scaleData(X,Y(:,k));
    [X_,Y_] = uniqueDataTol(X_,Y_,0.05);
%     X_ = unscaleData(X_,maxX);
%     Y_ = unscaleData(Y_,maxY);

    figure; hold on
    try
        Y_fit = net(X_);
        plot(Y_,Y_fit,'.')
    catch
        Y_fit = net(X_');
        plot(Y_,Y_fit,'.')
        plot([-1 1],[-1 1])

    end


    maxY(k) = maxY_;
end

n1 = models_{1};
n2 = models_{2};

models = @(x) [n1(x);n2(x)];

return


%
% [X,Y,maxX,maxY] = scaleData(X,Y);

% get unique data points
% tol = 0.008;
% [X,Y] = uniqueDataTol(X,Y,tol);

% size(X)
%
% Y1 = Y(:,1);
% Y2 = Y(:,2);



% return


% figure; hold on
% plot(data(iCase).time,data(iCase).states(:,1))
% plot(data(iCase).time,data(iCase).state_derivatives(:,1))

%close all

% iCase = 4;
% X = [data(iCase).states,data(iCase).inputs];
% Y = [data(iCase).state_derivatives];

% iY = 1;
% b = Y(:,iY)
% A = X;
% x = A\b;
%
% figure; hold on; plot(X*c); plot(Y(:,iY))

% return


%
% % [~,IA,~] = uniquetol([X,Y],0.01,'ByRows',true);
% %
% % Tu = data(iCase).time;
% % Tu = Tu(IA);
% % Xu = X(IA,:);
% % Yu = Y(IA,:);
% %
% % dataset1 = [Xu,Yu(:,1)];
% % dataset2 = [Xu,Yu(:,2)];
%
% % size(IA)
%
% % [W,IA,IC] = uniquetol(X(:,end-2),0.05);
%
% W = X(:,end-2);
% % [~,W_centers] = hist(W,25)
% W_centers = linspace(min(W),max(W),50);
% W_midpoint = W_centers(1:end-1)+diff(W_centers)/2;
%
% nw = length(W_centers)-1;
%
%
%
% figure; hold on
%
%
% [X,Y,maxX,maxY] = scaleData(X,Y);
%
%
% loadflag = 1;
%
% if loadflag
%
% mod = load('trained_nets.mat');
%     models = mod.models;
%
% else
%
% for k = 1:nw
%
%     %W_centers(k)
%
%     IC = (W_centers(k) <= W) & (W_centers(k+1) >= W);
%
%     Xw = X(IC,:);
%     Yw = Y(IC,:);
%
%     model = nntrain(Xw,Yw);
%
%     models{k} = model;
%
%     subplot(ceil(sqrt(nw)),ceil(sqrt(nw)),k); hold on
%     xm = model(Xw')';
%     ym = Yw;
%     plot(xm,ym,'.')
%     plot([-1 1],[-1 1],'color',[0.5 0.5 0.5])
%     xlim([min(xm,[],"all") max(xm,[],"all")])
%     ylim([min(ym,[],"all") max(ym,[],"all")])
%
%     text(0,0,string(W_centers(k)))
%
%
% end
%
% % save
% save('trained_nets.mat','models')
%
% end
% % plot histogram
% figure
% hist(W,W_centers)

%return

%% run a simulation test
%iCase = 6;

% case to simulate
iCase = 1;


W_centers = [];
% maxX = 1;
% maxY = 1;
% models = net;

tf = 110;
t_ind = (data(iCase).time < tf);
TSPAN = [data(iCase).time(1) tf];
Y0 = data(iCase).states(1,:);
OPTIONS = odeset('RelTol',1e-5,'AbsTol',1e-5);
U = griddedInterpolant(data(iCase).time,data(iCase).inputs,'spline');
[TOUT,YOUT] = ode45(@(t,x) deriv(t,x,U,models,maxX,maxY,W_centers),TSPAN,Y0,OPTIONS);

% plot
figure
subplot(2,1,1); hold on
plot(TOUT,YOUT(:,1),'r')
plot(data(iCase).time(t_ind),data(iCase).states(t_ind,1),'k')
title('Ptfm Pitch [deg]')


subplot(2,1,2); hold on
plot(TOUT,YOUT(:,2),'r')
plot(data(iCase).time(t_ind),data(iCase).states(t_ind,2),'k')
title('Gen Speed [rpm]')



function [net,tr,xi,ai,t] = nndynamic(X,Y)

X = tonndata(X,false,false);
T = tonndata(Y,false,false);

% Choose a Training Function
% For a list of all training functions type: help nntrain
% 'trainlm' is usually fastest.
% 'trainbr' takes longer but may be better for challenging problems.
% 'trainscg' uses less memory. Suitable in low memory situations.
trainFcn = 'trainlm';  % Levenberg-Marquardt backpropagation.

% Create a Time Delay Network
inputDelays = 1:2;
hiddenLayerSize = 10;
net = timedelaynet(inputDelays,hiddenLayerSize,trainFcn);

% Prepare the Data for Training and Simulation
% The function PREPARETS prepares timeseries data for a particular network,
% shifting time by the minimum amount to fill input states and layer
% states. Using PREPARETS allows you to keep your original time series data
% unchanged, while easily customizing it for networks with differing
% numbers of delays, with open loop or closed loop feedback modes.
[x,xi,ai,t] = preparets(net,X,T);

% Setup Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;

% Train the Network
[net,tr] = train(net,x,t,xi,ai);


end
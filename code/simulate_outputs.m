clc; clear; close all;

% get current directory

current_folder = pwd;

% add DanZs openfast output processing toolbox to path
tool_path = [current_folder,'/matlab-toolbox'];
addpath(genpath(tool_path))

% path to folder where the openfast outputs are
Out_path = 'data';

% prefix and suffix of output files
prefix = 'IEA_w_TMD_';
suffix = '.outb';

% get all the files in the directory
Out_files = dir(fullfile(Out_path,[prefix,'*',suffix]));

% get length
nLinCases = length(Out_files);

% numbering scheme bases on the number of files
if nLinCases <= 10
    numstring = '%01d';
else
    numstring = '%02d';
end

% reqired outputs
output_names = {'RtVAvgxh','GenTq','BldPitch1','PtfmPitch','GenSpeed','GenPwr'};
n_names = length(output_names);
sim_plot = 1;

filterflag = 1;

% go through each case
for iCase = 1:nLinCases
 
    switch suffix 
        case '.outb'
            [Channels, ChanName, ChanUnit, DescStr] = ReadFASTbinary(fullfile(Out_path,[prefix,num2str(iCase-1,'%01d'),suffix]));
        case '.mat'
            load([prefix,num2str(iCase-1,'%01d'),suffix]);
    end
    
    %go through each output name and plot trajectory
    if sim_plot
        figure(iCase)
        for idx = 1:n_names
            subplot(3,2,idx)
            ind = contains(ChanName,output_names{idx});
            %      find(ind)
            plot(Channels(:,1),Channels(:,ind))
            title([output_names{idx},'',ChanUnit{ind}])
            xlim([100,800])
        end
    end
    iTime = 1;
    iStates = [23,9]; %,22,24,25,26,27]; %
    iInputs = [2,268,6]; % RtVAvgxh, GenTq, BldPitch1

    % extract
    t = Channels(:,iTime);
    x = Channels(:,iStates);x(:,1) = rad2deg(x(:,1));
    u = Channels(:,iInputs);
    
    % filtering the wind input
    if filterflag
        t_f     = 1;
        dt      = t(2) - t(1);

        nb      = floor(t_f/dt);
        b       = ones(nb,1)/nb;

        u(:,1)     = filtfilt(b,1,u(:,1));
    end


    % assign
    data(iCase).time = t;
    data(iCase).states = x;
    data(iCase).inputs = u;

    % construct polynomial interpolates
    pp_x = spline(t,x');

    % integrate 
    pp_Ix = fnint(pp_x,ppval(pp_x,t(1)));

    % calculate integrated states
    data(iCase).states = x;

    % take derivatives of polynomial interpolates
    pp_Dx = fnder(pp_x,1);

    % calculate approximate derivatives
    Dx = ppval(pp_Dx,t);

    % assign derivatives
    data(iCase).state_derivatives = Dx';

end



% combine DLCs
Time = vertcat(data(:).time);
States = vertcat(data(:).states);
Inputs = vertcat(data(:).inputs);
State_derivatives = vertcat(data(:).state_derivatives);

% data for fitting
X = [States,Inputs];
Y = State_derivatives;

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



% [~,IA,~] = uniquetol([X,Y],0.01,'ByRows',true);
% 
% Tu = data(iCase).time;
% Tu = Tu(IA);
% Xu = X(IA,:);
% Yu = Y(IA,:);
% 
% dataset1 = [Xu,Yu(:,1)];
% dataset2 = [Xu,Yu(:,2)];

% size(IA)

% [W,IA,IC] = uniquetol(X(:,end-2),0.05);

W = X(:,end-2);
% [~,W_centers] = hist(W,25)
W_centers = linspace(min(W),max(W),50);
W_midpoint = W_centers(1:end-1)+diff(W_centers)/2;

nw = length(W_centers)-1;



figure; hold on 


[X,Y,maxX,maxY] = scaleData(X,Y);


loadflag = 1;

if loadflag
    
mod = load('trained_nets.mat');
    models = mod.models;
    
else
    
for k = 1:nw

    %W_centers(k)

    IC = (W_centers(k) <= W) & (W_centers(k+1) >= W);

    Xw = X(IC,:);
    Yw = Y(IC,:);

    model = nntrain(Xw,Yw);

    models{k} = model;

    subplot(ceil(sqrt(nw)),ceil(sqrt(nw)),k); hold on
    xm = model(Xw')';
    ym = Yw;
    plot(xm,ym,'.')
    plot([-1 1],[-1 1],'color',[0.5 0.5 0.5])
    xlim([min(xm,[],"all") max(xm,[],"all")])
    ylim([min(ym,[],"all") max(ym,[],"all")])

    text(0,0,string(W_centers(k)))


end

% save
save('trained_nets.mat','models')

end
% plot histogram 
figure
hist(W,W_centers)

%return

%% run a simulation test
%iCase = 6;

% case to simulate
iCase = 6;


tf = 800;
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


return

function dx = deriv(t,x,U,net,maxX,maxY,W_centers)

% get wind speed
u = U(t);

% get wind speed
w = u(1);

% determine which nueral net must be used
net_num = 1:length(W_centers);
n = ceil(interp1(W_centers,net_num,w));

% get sclaed input values
X = [x',u]./maxX;

% use net
net = net{n};

% evaluate derivative value from neural network
dx = net(X');

% unscale
dx = dx.*maxY';

end


%
function [X,Y] = uniqueDataTol(X,Y,tol)

% todo: scale data

% determine unique values
[~,IA,~] = uniquetol([X,Y],tol,'ByRows',true);

% extract unique indices
X = X(IA,:);
Y = Y(IA,:);

end

% 
function [X,Y,maxX,maxY] = scaleData(X,Y)

% maximum column values
maxX = max(abs(X),[],1);
maxY = max(abs(Y),[],1);

% avoid zeros
maxX(maxX==0) = 1;
maxY(maxY==0) = 1;

% scale
X = X./maxX;
Y = Y./maxY;

end

function X = unscaleData(X,maxX)

X = X.*maxX;

end

%
function net = nntrain(X,Y)

% 
x = X';
t = Y';

% Choose a Training Function
% For a list of all training functions type: help nntrain
% 'trainlm' is usually fastest.
% 'trainbr' takes longer but may be better for challenging problems.
% 'trainscg' uses less memory. Suitable in low memory situations.
trainFcn = 'trainbr';  % Levenberg-Marquardt backpropagation.

% Create a Fitting Network
hiddenLayerSize = 25;
net = fitnet(hiddenLayerSize,trainFcn);

net.trainParam.max_fail = 10000;
net.trainParam.mu_max = 1e20;
net.trainParam.min_grad = 1e-25;
net.trainParam.epochs = 200;

% Setup Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;

% Train the Network
[net,tr] = train(net,x,t,'UseParallel','yes');

% Test the Network
y = net(x);
e = gsubtract(t,y);
performance = perform(net,t,y);

% View the Network
% view(net)


end
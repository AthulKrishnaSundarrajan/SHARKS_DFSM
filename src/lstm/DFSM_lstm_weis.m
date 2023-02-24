clc;clear;close all;

root_path = which('INSTALL_DFSM'); % obtain full function path
data_path = fullfile(fileparts(root_path), 'data', filesep);
fulldata_path = fullfile(data_path,'transition_case'); % new data

%% Options for different simulation results
% prefix and suffix of output files
prefix = 'lin_'; % simulations without TTDspFA
% prefix = 'lin_'; % simulations with TTDspFA

suffix = '.outb';

%%
% get all the files in the directory
Out_files = dir(fullfile(fulldata_path,[prefix,'*',suffix]));

% get length
nLinCases = length(Out_files);
 %nLinCases = 1;

% numbering scheme bases on the number of files
if nLinCases <= 10
    numstring = '%01d';
else
    numstring = '%02d';
end

% required outputs
output_names = {'RtVAvgxh','GenTq','BldPitch1','PtfmPitch','GenSpeed'};
% output_names = {'RtVAvgxh','GenTq','BldPitch1','PtfmPitch','GenSpeed','TTDspFA'};

%req_states = {'PtfmPitch','GenSpeed'};
req_states = {'PtfmPitch','GenSpeed'};

req_controls = {'RtVAvgxh','GenTq','BldPitch1'};
n_names = length(output_names);iTime = 1;
sim_plot = 0;

wind = cell(nLinCases,1);

reduce_samples = ~true;

% go through each case
for iCase = 1:nLinCases

    switch suffix
        case '.outb'
            [Channels, ChanName, ChanUnit, DescStr] = ReadFASTbinary(fullfile(fulldata_path,[prefix,num2str(iCase-1,numstring),suffix]));
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
            xlim([Channels(1,1),Channels(end,1)])
        end
    end

    % fix to bug
    for i = 1:length(req_states)
        iStates(i) = find(contains(ChanName,req_states{i}));
    end

    % fix to bug
    for i = 1:length(req_controls)
        iInputs(i) = find(contains(ChanName,req_controls{i}));
    end

    % extract
    t = Channels(:,iTime);
    x = Channels(:,iStates); %x(:,1) = deg2rad(x(:,1));
    u = Channels(:,iInputs);

    if reduce_samples

        nt = 10000;
        t_ = linspace(t(1),t(end),nt);
        x = interp1(t,x,t_,'pchip');
        u = interp1(t,u,t_,'pchip');
        t = t_;

    end

    % assign
    data(iCase).time = t -min(t);
    data(iCase).states = x;
    data(iCase).inputs = u;
    data(iCase).state_names = ChanName(iStates);
    data(iCase).input_names = ChanName(iInputs);

    wind{iCase} = u(:,1);

   
end



% specify which states derivative is included
dindex = 1:length(req_states);

% approximate state derivatives
data = approximateStateDerivatives(data,dindex,1,3);

% names
state_names = req_states;
for k = 1:length(dindex)
    state_names{end+1} = ['d',req_states{dindex(k)}];
end

states = data(1).states;
controls = data(1).inputs;
%controls(end,:) = []; states(end,:) = [];

Input = [controls states]';
state_derivatives = data(1).state_derivatives;
Output = state_derivatives';

 time = data(1).time;
% time = time - min(time);
%time(end) = [];


% plot
hf = figure;
hf.Color = 'w';

subplot(5,1,1)
plot(time,controls(:,1))
title('Wind Speed')

subplot(5,1,2)
plot(time,controls(:,2))
title('GenTq')

subplot(5,1,3)
plot(time,controls(:,3))
title('BldPitch')

subplot(5,1,4)
plot(time,states(:,1))
title('PtfmPitch')

subplot(5,1,5)
plot(time,states(:,2))
title('GenSpeed')


% split data into training and testing

nt = length(time);

if length(data) >1
    
    nsplit = 100;

    indTrain = 1:nt; ntrail = length(indTrain);
    
    time_train = time;

    XTrainA = Input;
    TTrainA = Output;

    states2 = data(2).states;
    controls2 = data(2).inputs;
    
    XTestA = [controls2 states2]';
    state_derivatives2 = data(2).state_derivatives;
    TTestA = state_derivatives2';
    
    time_test = data(2).time;
    time_test = time_test - min(time);
    ntest = length(time_test);


else
    nsplit = 50; %<---------------------------

    indTrain = 1: floor(nsplit/100*nt);ntrain = length(indTrain);
    indTest = indTrain(end)+1:nt; ntest = length(indTest);
    
    XTrainA = Input(:,indTrain);
    TTrainA = Output(:,indTrain);
    
    XTestA = Input(:,indTest);
    TTestA = Output(:,indTest);
    
    time_train = time(indTrain); 
    time_test = time(indTest);
end
%------------------

hf = figure;
hf.Color = 'w';

ind = 1;

subplot(2,2,ind)
hold on;
plot(time_test,TTestA(ind,:),'k')
%plot(time_test,dx_lstm(ind,:),'r--')
xlabel('Time'); ylabel('dx_1')
%legend('original','DFSM LSTM')


ind = ind + 1;

subplot(2,2,ind)
hold on;
plot(time_test,TTestA(ind,:),'k')
%plot(time_test,dx_lstm(ind,:),'r--')
xlabel('Time'); ylabel('dx_2')
%legend('original','DFSM LSTM')


ind = ind + 1;

subplot(2,2,ind)
hold on;
plot(time_test,TTestA(ind,:),'k')
%plot(time_test,dx_lstm(ind,:),'r--')
xlabel('Time'); ylabel('dx_3')
%legend('original','DFSM LSTM')


ind = ind + 1;

subplot(2,2,ind)
hold on;
plot(time_test,TTestA(ind,:),'k')
%plot(time_test,dx_lstm(ind,:),'r--')
xlabel('Time'); ylabel('dx_4')
%legend('original','DFSM LSTM')

%-------------------------------------

n_channels = size(Input,1);
n_outputs = size(Output,1);

muX = mean(XTrainA,2);
sigmaX = std(XTrainA,0,2);

muT = mean(TTrainA,2);
sigmaT = std(TTrainA,0,2);

% reshape into cell array
ncells = 100;
indTrain_ = reshape(indTrain,[],ncells);


XTrain = cell(1,ncells);
TTrain = cell(1,ncells);

for i = 1:ncells
    ind_ = indTrain_(:,i);

    XTrain{i} = (XTrainA(:,ind_) - muX)./sigmaX;
    TTrain{i} = (TTrainA(:,ind_) - muT)./sigmaT;

end


% define architecture
layers = [
    sequenceInputLayer(n_channels)
    lstmLayer(128)
    fullyConnectedLayer(n_outputs)
    regressionLayer];

% training options
options = trainingOptions("adam", ...
    MaxEpochs=200, ...
    SequencePaddingDirection="left", ...
    Shuffle="every-epoch", ...
    Plots="training-progress", ...
    Verbose=1);

% train

tic;
net = trainNetwork(XTrain,TTrain,layers,options);
t_train = toc;

% normalize test
XTest = (XTestA - muX)./sigmaX;
TTest = (TTestA - muT)./sigmaT;

% test



dx_lstm = zeros(n_outputs,ntest);

% predict using the trained net
tic;
for i = 1:ntest

    % extract
    input_ = XTest(:,i);

    % predict
    [net,dx_lstm(:,i)] = predictAndUpdateState(net,input_);


end
t_predict = toc;

hf = figure;
hf.Color = 'w';

ind = 1;

subplot(2,2,ind)
hold on;
plot(time_test,TTest(ind,:),'k')
plot(time_test,dx_lstm(ind,:),'r--')
xlabel('Time'); ylabel('dx_1')
legend('original','DFSM LSTM')


ind = ind + 1;

subplot(2,2,ind)
hold on;
plot(time_test,TTest(ind,:),'k')
plot(time_test,dx_lstm(ind,:),'r--')
xlabel('Time'); ylabel('dx_2')
legend('original','DFSM LSTM')


ind = ind + 1;

subplot(2,2,ind)
hold on;
plot(time_test,TTest(ind,:),'k')
plot(time_test,dx_lstm(ind,:),'r--')
xlabel('Time'); ylabel('dx_3')
legend('original','DFSM LSTM')


ind = ind + 1;

subplot(2,2,ind)
hold on;
plot(time_test,TTest(ind,:),'k')
plot(time_test,dx_lstm(ind,:),'r--')
xlabel('Time'); ylabel('dx_4')
legend('original','DFSM LSTM')

%% run simulations

dx_lstm = dx_lstm.*sigmaT;
dx_lstm = dx_lstm + muT;

% state derivativs for validation
dx_test = state_derivatives(indTest,:);


dx_fun = @(t) interp1(time_test,dx_test,t,'pchip');
dx_dfsm = @(t) interp1(time_test,dx_lstm',t,'pchip');

state_test = states(indTest,:);
X0 = state_test(1,:)';

tspan = [time_test(1),time_test(end)];
OPTIONS = odeset('reltol',1e-8);

[T_OF,X_OF] = ode45(@(t,x)dx_function(t,x,dx_fun),tspan,X0,OPTIONS);
[T_dfsm,X_dfsm] = ode45(@(t,x)dx_function(t,x,dx_dfsm),T_OF,X0,OPTIONS);

%% save results

if reduce_samples
    suffix = 'red.mat';
else
    suffix = '.mat';
end

matname = fullfile(data_path,['lstm_weis_results_',num2str(nsplit),suffix]);

save(matname,'time_test','dx_lstm','TTest','XTest','time_train','XTrain','TTrain','muX','muT','sigmaX',...
    'sigmaT','time','states','controls','state_derivatives','t_train','t_predict',"nsplit",'T_OF','X_OF','T_dfsm','X_dfsm');


return


function dx = dx_function(t,x,dx_fun)

    dx = dx_fun(t);
    dx = dx';

end


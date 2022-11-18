close all; clear; clc
rng(34879)

% nonlinear derivative function
f = @(t,y,u,m) [y(2);m*(1 - y(1)^2)*y(2) - y(1) + u(1)];

% parameters
m = 1;

% time horizon
t0 = 0;
tf = 5;

% control (input)
U = @(t,p) 5*sin(exp(0.08.*t).*t + p);

% initial states and phase array
n = 100;
Y0array = -5 + 10*rand(2,n);
p = pi/2*rand(n,1);

% go through each simulation
for k = 1:size(Y0array,2)

    % extract initial states
    Y0 = Y0array(:,k);

    % simulation options
    OPTIONS = odeset('RelTol',1e-12,'AbsTol',sqrt(eps));

    % run nonlinear simulation
    [T1,Y1] = ode45(@(t,y) f(t,y,U(t,p(k)),m),[t0 tf],Y0,OPTIONS);

    % construct polynomial interpolate
    PY = spline(T1',Y1');

    % derivative of polynomial interpolate (state derivatives)
    PDY = fnder(PY,1);

    % compute approx state derivatives on simulation time mesh
    PDY1 = ppval(PDY,T1);

    % compute original nonlinear derivative values (for testing)
    DY1 = zeros(size(Y1'));
    for kk = 1:length(T1)
        DY1(:,kk) = f(T1(kk),Y1(kk,:),U(T1(kk),p(k)),m);

    end

    % assign
    T_array{k} = T1(:)';
  
    Iarray{k} = [Y1 U(T1,p(k))]';
    PDY1array{k} = PDY1;
    DY1array{k} = DY1;

end

% split into training and 
nTrain = 1:90;
nTest = 90+1:100;

% seperate training and test data sets
T_Train = T_array(nTrain);
XTrain = Iarray(nTrain);
TTrain = PDY1array(nTrain);

T_test = T_array(nTest);
XTest = Iarray(nTest);
TTest = PDY1array(nTest);

n_channels = size(Iarray{1},1);

% normalize
muX = mean(cat(2,XTrain{:}),2);
sigmaX = std(cat(2,XTrain{:}),0,2);

muT = mean(cat(2,TTrain{:}),2);
sigmaT = std(cat(2,TTrain{:}),0,2);

for n = 1:numel(XTrain)
    XTrain{n} = (XTrain{n} - muX) ./ sigmaX;
    TTrain{n} = (TTrain{n} - muT) ./ sigmaT;
end


% define architecture
layers = [
    sequenceInputLayer(n_channels)
    lstmLayer(128)
    fullyConnectedLayer(2)
    regressionLayer];

% training options
options = trainingOptions("adam", ...
    MaxEpochs=400, ...
    SequencePaddingDirection="left", ...
    Shuffle="every-epoch", ...
    Plots="training-progress", ...
    Verbose=1);

% train
net = trainNetwork(XTrain,TTrain,layers,options);

% normalize
for n = 1:5
    X = XTest{n};
    T = TTest{n};
    XTest{n} = (X - muX) ./ sigmaX;
    TTest{n} = (T - muT) ./ sigmaT;
end

%------------------------------------------------
% predict one test
idx = 2;


time = T_test{idx};
X = XTest{idx};
T = TTest{idx};

% initialize
PDY = zeros(2,length(time));

% predict values
for i = 1:length(time)

    xt = X(:,i);

    [net,PDY(:,i)] = predictAndUpdateState(net,xt);


end

% plot
hf = figure;
hf.Color = 'w';

ind = 1;
subplot(2,1,ind)
hold on;
plot(time,T(ind,:),'k')
plot(time,PDY(ind,:),'r')
xlabel('Time'); ylabel('dx_1')
legend('original','DFSM LSTM')


ind = 2;
subplot(2,1,ind)
hold on;
plot(time,T(ind,:),'k')
plot(time,PDY(ind,:),'r')
xlabel('Time'); ylabel('dx_2')
legend('original','DFSM LSTM')



%--------------------

% predict a second test
idx = 5;

% extract
time = T_test{idx};
X = XTest{idx};
T = TTest{idx};

% initialize
PDY = zeros(2,length(time));

% loop through and predict states
for i = 1:length(time)
    
    % extract
    xt = X(:,i);
    
    % predict
    [net,PDY(:,i)] = predictAndUpdateState(net,xt);


end

% plot

hf = figure;
hf.Color = 'w';

ind = 1;
subplot(2,1,ind)
hold on;
plot(time,T(ind,:),'k')
plot(time,PDY(ind,:),'r')
xlabel('Time'); ylabel('dx_1')
legend('original','DFSM LSTM')


ind = 2;
subplot(2,1,ind)
hold on;
plot(time,T(ind,:),'k')
plot(time,PDY(ind,:),'r')
xlabel('Time'); ylabel('dx_2')
legend('original','DFSM LSTM')


return


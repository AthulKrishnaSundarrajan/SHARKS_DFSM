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
n = 20;
Y0array = -5 + 10*rand(2,n);
p = pi/2*rand(n,1);

% go through each simulation
for k = 1:size(Y0array,2)

    % extract initial states
    Y0 = Y0array(:,k);

    % simulation options
    OPTIONS = odeset('RelTol',1e-6,'AbsTol',sqrt(eps));

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
    T1array{k} = T1(:)';
    Y1array{k} = Y1';
    U1array{k} = U(T1,p(k))';
    PDY1array{k} = PDY1;
    DY1array{k} = DY1;

end

% combine
T1 = horzcat(T1array{:});
Y1 = horzcat(Y1array{:});
U1 = horzcat(U1array{:});
PDY1 = horzcat(PDY1array{:});
DY1 = horzcat(DY1array{:});

% obtain unique inputs/outputs
[~,IA] = uniquetol([Y1;U1;PDY1]',1e-1,'ByRows',true);

% inputs
I = [Y1;U1]';
I = I(IA,:)';

% outputs
O = PDY1(:,IA);

%% phase plot of the two states
hf = figure; hold on
hf.Color = 'w';

% state 1 vs. state 2
plot(Y1(1,:),Y1(2,:),'.')
xlabel('Y1'); ylabel('Y2');

%% state-control plot
hf = figure; hold on
hf.Color = 'w';

% state 1 vs. state 2
plot3(I(1,:),I(2,:),I(3,:),'.')
xlabel('Y1'); ylabel('Y2'); zlabel('U');
title('Data Points for Fitting');

%% plot derivative
hf = figure; hold on
hf.Color = 'w';

% actual vs. approximate
plot(PDY1(:)-DY1(:),'.')
ylabel('error in DY');

% plot(T1(:),DY1','.')
% plot(T1(:),PDY1','.')

%% train the NN
% inputs
I = [Y1;U1];

% outputs
O = PDY1;

% choose a Training Function
% For a list of all training functions type: help nntrain
% 'trainlm' is usually fastest.
% 'trainbr' takes longer but may be better for challenging problems.
% 'trainscg' uses less memory. Suitable in low memory situations.
trainFcn = 'trainlm';  % Levenberg-Marquardt backpropagation.

% create a Fitting Network
hiddenLayerSize = 15;
net = fitnet(hiddenLayerSize,trainFcn);

% setup Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;

% train the Network
%[net,tr] = train(net,I,O);

% generate matlab function
%genFunction(net,'DFSM_TEST_net');

%% states
% initialize plot
hf = figure; hold on
hf.Color = 'w';

% go through each simulation
for k = 1:length(T1array)

    % initial states
    Y0 = Y0array(:,k);

    % run NN simulation
    % [T2,Y2] = ode45(@(t,y) net([y' U(t,p(k))]'),[t0 tf],Y0,OPTIONS);
    [T2,Y2] = ode45(@(t,y) DFSM_TEST_net([y' U(t,p(k))]'),T1array{k},Y0,OPTIONS);

    % initialize subplot
    subplot(length(T1array),1,k); hold on

    % plot the states
    plot(T1array{k},Y1array{k})
    plot(T2,Y2)

end

%% new initial point and control function
% initial states
Y02 = [-1;2];

% control (input)
U = @(t) exp(-t);

% run nonlinear simulation
[T1,Y1] = ode45(@(t,y) f(t,y,U(t),m),[t0 tf],Y02,OPTIONS);

% run NN simulation
[T2,Y2] = ode45(@(t,y) DFSM_TEST_net([y' U(t)]'),T1,Y02,OPTIONS);

% initialize plot
hf = figure; hold on
hf.Color = 'w';

% plot the states
subplot(3,1,1); hold on
plot(T1,U(T1));
xlabel('Time [s]')
ylabel('Controls')

subplot(3,1,2); hold on
plot(T1,Y1(:,1))
plot(T2,Y2(:,1),'-o')
legend('Original','DFSM')

xlabel('Time [s]')
ylabel('State 1')

subplot(3,1,3); hold on
plot(T1,Y1(:,2))
plot(T2,Y2(:,2),'-o')
legend('Original','DFSM')

xlabel('Time [s]')
ylabel('State 2')

%% new initial point and control function
% initial states
Y02 = [-5;5];

% control (input)
U = @(t) round(sin(t)) + exp(-0.1*t).*round(cos(t)) + 0.5;

% run nonlinear simulation
[T1,Y1] = ode45(@(t,y) f(t,y,U(t),m),[t0 tf],Y02,OPTIONS);

% run NN simulation
[T2,Y2] = ode45(@(t,y) DFSM_TEST_net([y' U(t)]'),[t0 tf],Y02,OPTIONS);

Y2_downsampled = @(t) interp1(T2,Y2,t,'pchip');
n = 100;
T1len = length(T1);

T2_ = linspace(T2(1),T2(end),n);

Y2_ = Y2_downsampled(T2_);

% initialize plot
close all;
hf = figure; hold on
hf.Color = 'w';
commonFigureProperties

% plot the states
%subplot(3,1,1); hold on
plot(T1,U(T1),'linewidth',1.5);
xlabel('Time [s]')
ylabel('u')

% axes
taskflag = 'axes'; commonFigureTasks;

% initialize plot
hf = figure; hold on
hf.Color = 'w';

plot(T1,Y1(:,1),'linewidth',1.5)
plot(T2_,Y2_(:,1),'o','linewidth',1.5)
legend('Original','DFSM')

xlabel('Time [s]')
ylabel('$x_1$')

% axes
taskflag = 'axes'; commonFigureTasks;

fontlegend = 15; nCol = 2;lcn = 'best';

taskflag = 'legend';commonFigureTasks;

% initialize plot
hf = figure; hold on
hf.Color = 'w';

plot(T1,Y1(:,2),'linewidth',1.5)
plot(T2_,Y2_(:,2),'o','linewidth',1.5)
legend('Original','DFSM')

xlabel('Time [s]')
ylabel('$x_2$')

% axes
taskflag = 'axes'; commonFigureTasks;

fontlegend = 15; nCol = 2;lcn = 'best';

taskflag = 'legend';commonFigureTasks;
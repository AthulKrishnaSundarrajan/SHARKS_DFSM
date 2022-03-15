function [mdl,X,Y] = createDFSM(data,linearModels)

test = 'LTI';

% combine DLCs
Time = vertcat(data(:).time);
States = vertcat(data(:).states);
Inputs = vertcat(data(:).inputs);
State_derivatives = vertcat(data(:).state_derivatives);

% data for fitting
X = [States,Inputs];
Y = State_derivatives;

%
iCase = 1;
U = griddedInterpolant(data(iCase).time,data(iCase).inputs,'spline');
dGenSpeed = griddedInterpolant(data(iCase).time,data(iCase).state_derivatives(:,2),'spline');

%
t0 = 100;
td = 500;
TSPAN = [t0 t0+td];
Y0 = data(iCase).states(t0 == data(iCase).time,:);
OPTIONS = odeset('RelTol',1e-5,'AbsTol',1e-5);

%
switch test
    case 'LTI'
        [Ai,Bi,Di] = DFSM_LTI(Inputs,States,State_derivatives);
        [TOUT,YOUT] = ode45(@(t,x) deriv_single(t,x,U,Ai,Bi,Di),TSPAN,Y0,OPTIONS);
        mdl = [];

    case 'LPV'
        [Ai,Bi,Di] = DFSM_LPV(Inputs,States,State_derivatives);
        [TOUT,YOUT] = ode45(@(t,x) deriv1(t,x,U,Ai,Bi,Di,dGenSpeed),TSPAN,Y0,OPTIONS);

    case 'NN'
        %
        [X_,Y_,maxX,maxY] = scaleData(X,Y);
        [X_,Y_] = uniqueDataTol(X_,Y_,0.02);

        root_path = which('INSTALL_DFSM'); % obtain full function path
        data_path = fullfile(fileparts(root_path), 'data', filesep);
        fulldata_path = fullfile(data_path,'idetc');

        net = nntrain(X_,Y_(:,[2,4,5]));
        genFunction(net,fullfile(fulldata_path,'mynet.m'))

        [TOUT,YOUT] = ode45(@(t,x) deriv_net2(t,x,U,maxX,maxY,dGenSpeed),TSPAN,Y0,OPTIONS);

%
% [TOUT,YOUT] = ode45(@(t,x) deriv_idw(t,x,U,X_,Y_),TSPAN,Y0,OPTIONS);

% [TOUT,YOUT] = ode45(@(t,x) deriv_net(t,x,U,[],[],dGenSpeed),TSPAN,Y0,OPTIONS);

%

end

% figure
% try
%     Y_fit = net(X_);
%     plot(Y_,Y_fit,'.')
% catch
%     Y_fit = net(X_')';
%     plot(Y_,Y_fit,'.')
% end
%
%
% mdl = [];
% return


% scale data

% offset = 0.2;

% plot
hf = figure;
hf.Color = 'w';
hf.Position = [1000 700 560 640];

subplot(5,1,1); hold on
plot(TOUT,YOUT(:,1),'r')
plot(data(iCase).time,data(iCase).states(:,1),'k')
title('Ptfm Pitch')
xlim([t0 t0+td])

subplot(5,1,2); hold on
plot(TOUT,YOUT(:,2),'r')
plot(data(iCase).time,data(iCase).states(:,2),'k')
title('Gen Speed')
xlim([t0 t0+td])

subplot(5,1,3); hold on
plot(TOUT,YOUT(:,3),'r')
plot(data(iCase).time,data(iCase).states(:,3),'k')
title('TTDspFA')
xlim([t0 t0+td])

subplot(5,1,4); hold on
plot(TOUT,YOUT(:,4),'r')
plot(data(iCase).time,data(iCase).states(:,4),'k')
title('Deriv Ptfm Pitch')
xlim([t0 t0+td])

subplot(5,1,5); hold on
plot(TOUT,YOUT(:,5),'r')
plot(data(iCase).time,data(iCase).states(:,5),'k')
title('Deriv TTDspFA')
xlim([t0 t0+td])

end

function dx = deriv_idw(t,x,U,X,Y)

% get wind speed
u = U(t);

xu = [x;u(:)];

p = 4;
rad = 100;
L = 1;

dx(1) = x(4);
dx(2) = idw(X,Y(:,2),xu',p,rad,L);
dx(3) = x(5);
dx(4) = idw(X,Y(:,4),xu',p,rad,L);
dx(5) = idw(X,Y(:,5),xu',p,rad,L);

disp(t)

dx = dx(:);

end

function dx = deriv_net(t,x,U,maxX,maxY,dGenSpeed)

% get wind speed
u = U(t);

xu = [x(:);u(:)];

%
dx = mynet(xu);

dx(1) = x(4);
dx(3) = x(5);
% dx(2) = dGenSpeed(t);

disp(t)

end

function dx = deriv_net2(t,x,U,maxX,maxY,dGenSpeed)

% get wind speed
u = U(t);

xu = [x(:);u(:)];
xu = xu./maxX(:);

maxY = maxY([2,4,5]);

%
dx_net = mynet(xu);
dx_net = dx_net.*maxY(:);

dx = zeros(5,1);

dx(1) = x(4);
dx(2) = dx_net(1);
dx(3) = x(5);
dx(4) = dx_net(2);
dx(5) = dx_net(3);
% dx(2) = dGenSpeed(t);

disp(t)

end

function dx = deriv1(t,x,U,Ai,Bi,Di,dGenSpeed)

% get wind speed
u = U(t);
w = u(1);
u = u(1:3);

if w > 16
    w = 16;
end

%
dx = Ai(w)*x + Bi(w)*u(:) + Di(w);

dx(1) = x(4);
dx(3) = x(5);

% dx(2) = dGenSpeed(t);

end

function dx = deriv_single(t,x,U,Ai,Bi,Di)

% get wind speed
u = U(t);
w = u(1);
% u = u;

dx = Ai*x + Bi*u(:) + Di;

dx(1) = x(4);
dx(3) = x(5);

% disp(t)

end

function [Ai,Bi,Di] = DFSM_LTI(Inputs,States,State_derivatives)

% construst data inputs
data = [States,Inputs,zeros(size(States,1),1)];

% find the best linear model
ABD = linsolve(data,State_derivatives);

% transpose
ABD = ABD';

% extract
Ai = ABD(:,1:5);
Bi = ABD(:,6:8);
Di = ABD(:,9);

end


function [Ai,Bi,Di] = DFSM_LPV(Inputs,States,State_derivatives)

w = Inputs(:,1);

W = linspace(7,16,50);
dw = W(2)-W(1);
offset = 3;

% W = 12;
% offset = 1;

% W = mean(w);
% offset = 100;

maxW = max(w);
minW = min(w);

for k = 1:length(W)

    m = W(k);

    if m+offset > maxW
        offset = maxW - m;
    end

    if m-offset < minW
        offset = m-minW;
    end

    % find all values within small range
    Iw = (m+offset >= w) & (m-offset <= w);

    % extract
    x_ = States(Iw,:);
    u_ = Inputs(Iw,:);
%     u_(:,1) = [];

    Dx_ = State_derivatives(Iw,:);

    in = [x_,u_,zeros(sum(Iw),1)];

    AB = linsolve(in,Dx_);

    mdl_ = fitlm([x_,u_],Dx_(:,2));

    Beta = linsolve([x_,u_],Dx_);

    AB = AB';
%     mdl(k).A = [0 0 0 1 0; AB(2,1:5); 0 0 0 0 1; AB(1,1:5); AB(3,1:5);];
%     mdl(k).B = [0 0 0; AB(2,6:8); 0 0 0; AB(1,6:8); AB(3,6:8)];
%     mdl(k).D = [0; AB(2,9); 0; AB(1,9); AB(3,9)];

    mdl(k).A = AB(:,1:5);
    mdl(k).B = AB(:,6:8);
    mdl(k).D = AB(:,9);


%     {'PtfmPitch','GenSpeed','TTDspFA','d_PtfmPitch','d_TTDspFA'};
%     {'d_PtfmPitch','d_GenSpeed','d_TTDspFA','d2_PtfmPitch','d2_TTDspFA'};
%     {'d2_PtfmPitch','d_GenSpeed','d2_TTDspFA'}


    Dx_mdl = [x_,u_,ones(sum(Iw),1)]*AB';

    mdl(k).mse = sum((Dx_-Dx_mdl).^2,'all')/numel(Dx_);

%     figure; hold on
%     plot(Dx_,Dx_mdl,'.')
%     axis manual
%     plot([-1 1],[-1 1],'k')

    eig(mdl(k).A)

end

% mdl = [];


for k = 1:length(mdl)

    B(:,k) = mdl(k).B(:);
    A(:,k) = mdl(k).A(:);

end

figure; plot(W,A)
figure; plot(W,B)

for k = 1:length(mdl)
    A_interp(k,:,:) = mdl(k).A;
    B_interp(k,:,:) = mdl(k).B;
    D_interp(k,:,:) = mdl(k).D;
end

%
A_pp = interp1(W,A_interp,'nearest','pp');
Ai = @(w) ppval(A_pp,w);
B_pp = interp1(W,B_interp,'nearest','pp');
Bi = @(w) ppval(B_pp,w);
D_pp = interp1(W,D_interp,'nearest','pp');
Di = @(w) ppval(D_pp,w);

% Ai = mdl(1).A;
% Bi = mdl(1).B;
% Di = mdl(1).D;


end
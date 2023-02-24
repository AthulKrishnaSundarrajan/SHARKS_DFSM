function model = createDFSM(data,model,iCase,state_names,input_names)

% extract state and input information
Time = vertcat(data(iCase).time);
States = vertcat(data(iCase).states);
Inputs = vertcat(data(iCase).inputs);
State_derivatives = vertcat(data(iCase).state_derivatives);
exportflag = 0;
% data for fitting
X = [States,Inputs];
Y = State_derivatives;

% number of states
ns = size(States,2);

% input
U = griddedInterpolant(data(iCase).time,data(iCase).inputs,'spline');
% dGenSpeed = griddedInterpolant(data(iCase).time,data(iCase).state_derivatives(:,2),'spline');

% simulation options
t0 = min(Time);
td = max(Time)-min(Time);
TSPAN = [t0 t0+td];
Y0 = data(iCase).states(t0 == Time,:);
OPTIONS = odeset('RelTol',1e-5,'AbsTol',1e-5);

% extract options
sim_type = model.sim_type;
mdl = model.mdl;

mdlfalg = isempty(mdl);
% generate models

switch sim_type

    case 'LTI'

        % create or load model
        if isempty(mdl)
            [Ai,Bi,Di] = DFSM_LTI(Inputs,States,State_derivatives);
            model.mdl = {Ai,Bi,Di};
        else
            Ai = mdl{1}; Bi = mdl{2}; Di = mdl{3};
        end

        % run simulation
         [TOUT,YOUT] = ode45(@(t,x) deriv_lti(t,x,U,Ai,Bi,Di),TSPAN,Y0,OPTIONS);

    case 'LPV'

        wmin = min(Inputs(:,1))+1;
        wmax = max(Inputs(:,1))-1;

        % create or load model
        if isempty(mdl)
            [Ai,Bi,Di] = DFSM_LPV(Inputs,States,State_derivatives,wmin,wmax);
            model.mdl = {Ai,Bi,Di};
        else
            Ai = mdl{1}; Bi = mdl{2}; Di = mdl{3};
        end

        % run simulation
        [TOUT,YOUT] = ode45(@(t,x) deriv_lpv(t,x,U,Ai,Bi,Di,wmin,wmax),TSPAN,Y0,OPTIONS);

    case 'NN'

        % create or load model
        if isempty(mdl)
            [X_,Y_,maxX,maxY] = scaleData(X,Y);
            [X_,Y_,I] = uniqueDataTol(X_,Y_,0.00);
            X = X(I,:);
            Y = Y(I,:);

            root_path = which('INSTALL_DFSM'); % obtain full function path
            data_path = fullfile(fileparts(root_path), 'data', filesep);
            fulldata_path = fullfile(data_path,'DFSM_wo_FA2');

            net = nntrain(X,Y);
            genFunction(net,fullfile(fulldata_path,'mynet.m'))
            addpath(fulldata_path)
            model.mdl = {maxX,maxY};

        else
            maxX = model.mdl{1}; maxY = model.mdl{2};

        end

        % run simulation
        [TOUT,YOUT] = ode45(@(t,x) deriv_net2(t,x,U,maxX,maxY,[]),TSPAN,Y0,OPTIONS);
end

% plot
hf = figure;
hf.Color = 'w';
hf.Position = [1000 700 560 640];

if mdlfalg

    title_name = [sim_type,' - Training'];
    fig1_name = 'sim_training';
    fig2_name = 'power_training';

else
    title_name = [sim_type,' - Validation'];
    fig1_name = 'sim_validation';
    fig2_name = 'power_validation';
end

sgtitle(title_name)
nl = length(state_names);

for i = 1:nl

subplot(nl,1,i); hold on
plot(TOUT,YOUT(:,i),'r')
plot(data(iCase).time,data(iCase).states(:,i),'k')
title(state_names{i})
xlim([t0 t0+td])

end


if exportflag
    savename = fig1_name;
    pathpdf = mfoldername(mfilename('fullpath'),'final');
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end


% plot power
hf = figure; hold on
hf.Color = 'w';

% extract genspeed and gentorq
IGenSpeed = strcmp(state_names,'GenSpeed');
GenSpeed = YOUT(:,IGenSpeed);
UOUT = U(TOUT);
IGenTq = strcmp(input_names,'GenTq');
GenTq = UOUT(:,IGenTq);
Iwind = strcmp(input_names,'RtVAvgxh');


plot(TOUT,GenSpeed.*GenTq,'r')
plot(Time,States(:,IGenSpeed).*Inputs(:,IGenTq),'k')
xlabel('Time [s]')
ylabel('Generator Power [kW]')
legend('DFSM','OpenFAST')
%ylim([0 2e5])
sgtitle(title_name)

if exportflag
    savename = fig2_name;
    pathpdf = mfoldername(mfilename('fullpath'),'final');
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end






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
dx(2) = x(5);
dx(3) = x(6);

disp(t)

end

function dx = deriv_net2(t,x,U,maxX,maxY,dGenSpeed)

% get wind speed
u = U(t);

xu = [x(:);u(:)];
%xu = xu./maxX(:);

% maxY = maxY([2,4,5]);

%
dx_net = mynet(xu);
dx = dx_net;

% ns = length(dx_net);
% 
% dx = zeros(2*ns,1);
% 
% dx(1:ns) = x(ns+1:end);
% dx(ns+1:end) = dx_net;

end

function dx = deriv_lti(t,x,U,Ai,Bi,Di)

% get inputs
u = U(t);

% number of states
ns = length(Ai)/2;

% compute state derivatives (LTI model)
dx = Ai*x + Bi*u(:)+ Di;

% update known derivative values
%dx(1:ns) = x(ns+1:end);

% disp(t)

end

function dx = deriv_lpv(t,x,U,Ai,Bi,Di,wmin,wmax)

% get wind speed
u = U(t);
w = u(1);
u = u(1:3);

if w > wmax
    w = wmax;
elseif w < wmin
    w = wmin;
end

%
dx = Ai(w)*x + Bi(w)*u(:); % + Di(w);

% number of states
ns = length(Ai(w))/2;

% update known derivative values
%dx(1:ns) = x(ns+1:end);

end

function [Ai,Bi,Di] = DFSM_LTI(Inputs,States,State_derivatives)

% construst data inputs
data = [States,Inputs,zeros(size(States,1),1)];
%data = [States,Inputs];
% data = data./mean(data);

% find the best linear model
ABD = linsolve(data,State_derivatives);

% transpose
ABD = ABD';

% number of states and inputs
ns = size(States,2);
nu = size(Inputs,2);

% extract
Ai = ABD(:,1:ns);%Ai(abs(Ai)<1e-10) = 0;
Bi = ABD(:,(ns+1):(ns+nu));
Di = ABD(:,end);

end

function [Ai,Bi,Di] = DFSM_LPV(Inputs,States,State_derivatives,wmin,wmax)

w = Inputs(:,1);
maxW = max(w);
minW = min(w);

W = linspace(wmin,wmax,50);
dw = W(2)-W(1);
offset = 3;
IW_ = zeros(length(W),1);
Aeig = zeros(length(W),1);
% W = 12;
% offset = 1;

% W = mean(w);
% offset = 100;

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
    IW_(k) = sum(Iw);

    % extract
    x_ = States(Iw,:);
    u_ = Inputs(Iw,:);
%     u_(:,1) = [];

    Dx_ = State_derivatives(Iw,:);

    %in = [x_,u_,zeros(sum(Iw),1)];
    in = [x_,u_];
    AB = linsolve(in,Dx_);

    %mdl_ = fitlm([x_,u_],Dx_(:,2));

    %Beta = linsolve([x_,u_],Dx_);

    AB = AB';
%     mdl(k).A = [0 0 0 1 0; AB(2,1:5); 0 0 0 0 1; AB(1,1:5); AB(3,1:5);];
%     mdl(k).B = [0 0 0; AB(2,6:8); 0 0 0; AB(1,6:8); AB(3,6:8)];
%     mdl(k).D = [0; AB(2,9); 0; AB(1,9); AB(3,9)];

    % number of states and inputs
    ns = size(States,2);
    nu = size(Inputs,2);

    mdl(k).A = AB(:,1:ns);
    mdl(k).B = AB(:,(ns+1):(ns+nu));
    mdl(k).D = [];%AB(:,end);


%     {'PtfmPitch','GenSpeed','TTDspFA','d_PtfmPitch','d_TTDspFA'};
%     {'d_PtfmPitch','d_GenSpeed','d_TTDspFA','d2_PtfmPitch','d2_TTDspFA'};
%     {'d2_PtfmPitch','d_GenSpeed','d2_TTDspFA'}

    Dx_mdl = [x_,u_]*AB';

    mdl(k).mse = sum((Dx_-Dx_mdl).^2,'all')/numel(Dx_);

%     figure; hold on
%     plot(Dx_,Dx_mdl,'.')
%     axis manual
%     plot([-1 1],[-1 1],'k')

   Aeig(k) = sum(real(eig(mdl(k).A)) > 0);


end

% mdl = [];


for k = 1:length(mdl)

    B(:,k) = mdl(k).B(:);
    A(:,k) = mdl(k).A(:);

end

% figure; plot(W,A)
% figure; plot(W,B)

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
%D_pp = interp1(W,D_interp,'nearest','pp');
Di = []; % @(w) ppval(D_pp,w);

% Ai = mdl(1).A;
% Bi = mdl(1).B;
% Di = mdl(1).D;


end
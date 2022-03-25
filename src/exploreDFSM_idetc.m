clc; clear; close all;

root_path = which('INSTALL_DFSM'); % obtain full function path
data_path = fullfile(fileparts(root_path), 'data', filesep);
fulldata_path = fullfile(data_path,'idetc');

% load the simulation data
type =  'turb'; % 'turb';

% data folder
outpath = 'data';

% load matrices
switch type
    case 'turb'
        data = load(fullfile(fulldata_path,'NL_Turb.mat'));
        ChanName = data.NL_Turb_ChanName;
        Channels = data.NL_Turb_Channels;
        ChanUnit = data.NL_Turb_ChanUnit;

    case 'step'
        data = load(fullfile(fulldata_path,'NL_Step.mat'));
        ChanName = data.NL_Step_ChanName;
        Channels = data.NL_Step_Channels;
        ChanUnit = data.NL_Step_ChanUnit;

end

% Time
Time = Channels(:,1);

% required states
req_states = {'PtfmPitch','GenSpeed','TTDspFA'};
req_states = {'PtfmPitch','TTDspFA','GenSpeed'};

state_ind = find(contains(ChanName,req_states));

for i = 1:length(req_states)
    state_ind(i) = find(contains(ChanName,req_states{i}));
end

States = Channels(:,state_ind);

% load specific states
PtfmPitch = States(:,1); TTDspFA = States(:,2); GenSpeed = States(:,3);

% required controls
req_controls = {'RtVAvgxh','GenTq','BldPitch1'};
controls_ind = find(contains(ChanName,req_controls));
Controls = Channels(:,controls_ind);

% load specific controls
WindSpeed = Controls(:,1); GenTq = Controls(:,2); BldPitch = Controls(:,3);

clear data
data(1).time = Time;
data(1).states = States;
data(1).inputs = Controls;

dindex = [1 3];


data = approximateStateDerivatives(data,dindex);

state_names = req_states; state_names{end+1} = ['d',req_states{dindex(1)}]; state_names{end+1} = ['d',req_states{dindex(2)}];
% load linearized model
LinearModels = load(fullfile(fulldata_path,'pd_1.0_linear.mat'));

iCase = 1;
model.sim_type = 'LPV';
model.mdl = {};

model = createDFSM(data,model,iCase,state_names);

% [X_,Y_,maxX,maxY] = scaleData(X,Y);
% [X_,Y_] = uniqueDataTol(X_,Y_,0.15);
% X_ = unscaleData(X_,maxX);
% Y_ = unscaleData(Y_,maxY);
% 
% Y2 = Y_(:,2);
% Y22 = Y(:,2);
return

% combine DLCs
Time = vertcat(data(:).time);
States = vertcat(data(:).states);
Inputs = vertcat(data(:).inputs);
State_derivatives = vertcat(data(:).state_derivatives);

% data for fitting
X = [States,Inputs];
Y = State_derivatives;

k = 1;

[X_,Y_,maxX,maxY_] = scaleData(X,Y(:,k));
[X_,Y_] = uniqueDataTol(X_,Y_,0.15);

% net = newrb(X_',Y_',0,10,400);
net = nntrain(X_,Y_);


figure; hold on
Y_fit = net(X_')';
plot(Y_,Y_fit,'.')
plot([-1 1],[-1 1])

figure; hold on
[X_,Y_,maxX,maxY_] = scaleData(X,Y(:,k));
[X_,Y_] = uniqueDataTol(X_,Y_,0.05);
Y_fit = net(X_')';
plot(Y_,Y_fit,'.')
plot([-1 1],[-1 1])

return

% load linearized model
LinearModels = load(fullfile(outpath,'pd_1.0_linear.mat'));

% Set length
nl = length(LinearModels.P);

% go through each linear model
for iSys = 1:nl

    % extract
    sys = LinearModels.P{iSys};

    % number of states
    nStates(iSys) = size(LinearModels.P{iSys}.A,1);

    % wind speed
    w_ops(iSys,1) = LinearModels.WindSpeed(iSys);

    % operating points
    x_opsM(:,iSys) = LinearModels.SS_Ops(iSys).xop;
    u_opsM(:,iSys) = LinearModels.SS_Ops(iSys).uop;

    % matrices
    Am(:,:,iSys) = LinearModels.P{iSys}.A;
    Bm(:,:,iSys) = LinearModels.P{iSys}.B;

    % state indices
    StateNames = sys.StateName;
    iGenSpeed = find(strcmpi('ED First time derivative of Variable speed generator DOF (internal DOF index = DOF_GeAz), rad/s',StateNames)); % generator speed
    iPltPitch = find(strcmpi('ED Platform pitch tilt rotation DOF (internal DOF index = DOF_P), rad',StateNames)); % platform pitch
    iTTDspFA = find(strcmpi('ED 1st tower fore-aft bending mode DOF (internal DOF index = DOF_TFA1), m',StateNames));


end



return
clc; clear; close all;

% load the simulation data
type =  'turb'; % 'turb';

% data folder
outpath = 'data';

% load matrices
switch type
    case 'turb'
        data = load(fullfile(outpath,'NL_Turb.mat'));
        ChanName = data.NL_Turb_ChanName;
        Channels = data.NL_Turb_Channels;
        ChanUnit = data.NL_Turb_ChanUnit;
        
    case 'step'
        data = load(fullfile(outpath,'NL_Step.mat'));
        ChanName = data.NL_Step_ChanName;
        Channels = data.NL_Step_Channels;
        ChanUnit = data.NL_Step_ChanUnit;
        
end

% Time
Time = Channels(:,1);

% required states
req_states = {'PtfmPitch','TTDspFA','GenSpeed'};
state_ind = find(contains(ChanName,req_states));
States = Channels(:,state_ind);

% load specific states
PtfmPitch = States(:,1);TTDspFA = States(:,2); GenSpeed = States(:,3);

% required controls
req_controls = {'RtVAvgxh','GenTq','BldPitch1'};
controls_ind = find(contains(ChanName,req_controls));
Controls = Channels(:,controls_ind);

% load specific controls
WindSpeed = Controls(:,1); GenTq = Controls(:,2); BldPitch = Controls(:,3);


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
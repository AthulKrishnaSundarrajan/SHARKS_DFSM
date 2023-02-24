clc; clear; close all;

root_path = which('INSTALL_DFSM'); % obtain full function path
data_path = fullfile(fileparts(root_path), 'data', filesep);
fulldata_path = fullfile(data_path,'idetc');

% load the simulation data
type1 =  'turb'; % 'turb';
req_states = {'PtfmPitch','TTDspFA','GenSpeed'};
req_controls = {'RtVAvgxh','GenTq','BldPitch1'};

data = [];
data = loadIDETCdata(fulldata_path,type1,req_states,req_controls);

% validation
% type2 = 'step';
% data2 = loadIDETCdata(fulldata_path,type2,req_states,req_controls);
% 
% data(2) = data2;
% 
 dindex = 1:length(req_states);


data = approximateStateDerivatives(data,dindex,0,5);


s2 = smoothData2(data.time,data.states(:,2));

hf = figure;
hf.Color = 'w';hold on
plot(data.time,data.states(:,2),'r--')
plot(data.time,s2,'ko')

state_names = req_states;
for k = 1:length(dindex)
    state_names{end+1} = ['d',req_states{dindex(k)}];
end
input_names = req_controls;

% load linearized model
LinearModels = load(fullfile(fulldata_path,'pd_1.0_linear.mat'));

iCase = 1;
model.sim_type = 'LPV';
model.mdl = {};

model = createDFSM(data,model,iCase,state_names,input_names);

iCase = 2;
model = createDFSM(data,model,iCase,state_names,req_controls);

controls = vertcat(data(:).inputs);

wind = controls(:,1); wind = reshape(wind,[length(wind)/2,2]);

hf = figure;
hf.Color = 'w';

hold on;
histogram(wind(:,1))
histogram(wind(:,2));
legend({[type1,' (T)'],[type2,' (V)']})


tcase_range = [max(wind(:,1)),min(wind(:,1))];

hf = figure;
hf.Color = 'w';

hold on;

yline(tcase_range);
plot(data(1).time,wind(:,1))
plot(data(2).time,wind(:,2));
xlabel('Time [s]')
ylabel('Wind Speed [m/s]')

legend({'Training range','','Training','Validation'},'location','best')


return


function data = loadIDETCdata(fulldata_path,type,req_states,req_controls)


% load matrices
switch type
    case 'turb'
        sim_results = load(fullfile(fulldata_path,'NL_Turb.mat'));
        ChanName = sim_results.NL_Turb_ChanName;
        Channels = sim_results.NL_Turb_Channels;
        ChanUnit = sim_results.NL_Turb_ChanUnit;

    case 'step'
        sim_results = load(fullfile(fulldata_path,'NL_Step.mat'));
        ChanName = sim_results.NL_Step_ChanName;
        Channels = sim_results.NL_Step_Channels;
        ChanUnit = sim_results.NL_Step_ChanUnit;

end

% Time
Time = Channels(:,1);

if strcmpi(type,'step')

    tind = Time <= 600;

    Time = Time(tind);
    Channels = Channels(tind,:);

end

% required states
% req_states = {'PtfmPitch','GenSpeed','TTDspFA'};


state_ind = find(contains(ChanName,req_states));

for i = 1:length(req_states)
    state_ind(i) = find(contains(ChanName,req_states{i}));
end

States = Channels(:,state_ind);

% load specific states
PtfmPitch = States(:,1); TTDspFA = States(:,2); GenSpeed = States(:,3);

% required controls

for i = 1:length(req_controls)
    iInputs(i) = find(contains(ChanName,req_controls{i}));
end
u = Channels(:,iInputs);


% % load specific controls
% WindSpeed = u(:,1); GenTq = u(:,2); BldPitch = u(:,3);

data.time = Time;
data.states = States;
data.inputs = u;




end

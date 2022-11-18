% script to run DFSM LTI with PCA

clc; clear; close all

root_path = which('INSTALL_DFSM'); % obtain full function path
data_path = fullfile(fileparts(root_path), 'data', filesep);
fulldata_path = fullfile(data_path,'sim_17ms_600'); % new data

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
req_states = {'PtfmPitch','TTDspFA','GenSpeed'};

req_controls = {'RtVAvgxh','GenTq','BldPitch1'};
n_names = length(output_names);iTime = 1;
sim_plot = 0;

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

    % assign
    data(iCase).time = t;
    data(iCase).states = x;
    data(iCase).inputs = u;
    data(iCase).state_names = ChanName(iStates);
    data(iCase).input_names = ChanName(iInputs);

   
end

% specify which states derivative is included
dindex = 1:length(req_states);

% approximate state derivatives
data = approximateStateDerivatives(data,dindex,1,0.5);

% perform PCA and scale
data = performPCA(data);

% add state names
LinearModels = [];
state_names = req_states;
for k = 1:length(dindex)
    state_names{end+1} = ['d',req_states{dindex(k)}];
end

% generate model using first simulation
iCase = 1;
model.sim_type = 'LTI';
model.mdl = {};


function data = performPCA(data)

% fitting data
ind = 1;

% extract staes and inputs
states = data(ind).states;
inputs = data(ind).inputs;

% combine and scale using maximum value
fitdata = [states,inputs];
scale_max = max(fitdata);
fitdata_scaled = fitdata./scale_max;

% add to data structure
data(ind).fitdata_scaled = fitdata_scaled;
data(ind).scale_max = scale_max;

% perform PCA
[coeff,score,latent,tsquared,explained,mu] = pca(fitdata_scaled);

% store coeff and mu
data(ind).pca_coeff = coeff;
data(ind).mu_pca = mu;


end
function data = loadData_IEA_w_TMD

close all;

root_path = which('INSTALL_DFSM'); % obtain full function path
data_path = fullfile(fileparts(root_path), 'data', filesep);
fulldata_path = fullfile(data_path,'IEA_w_TMD2');

% prefix and suffix of output files
prefix = 'lin_'; % near rated region
% prefix = 'lin_rated_'; % rated
suffix = '.outb';

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
output_names = {'RtVAvgxh','GenTq','BldPitch1','PtfmPitch','GenSpeed','GenPwr'};
req_states = {'PtfmPitch','TipDxb1','GenSpeed'};
req_controls = {'RtVAvgxh','GenTq','BldPitch1'};
n_names = length(output_names);iTime = 1;
sim_plot = 1;

% go through each case
for iCase = 1:nLinCases

    switch suffix
        case '.outb'
            [Channels, ChanName, ChanUnit, DescStr] = ReadFASTbinary(fullfile(fulldata_path,[prefix,num2str(iCase-1,'%01d'),suffix]));
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
dindex = [1,3];

% approximate state derivatives
data = approximateStateDerivatives(data,dindex);

% add state names
LinearModels = [];
state_names = req_states; state_names{end+1} = ['d',req_states{dindex(1)}]; state_names{end+1} = ['d',req_states{dindex(2)}];

% generate model using first simulation
iCase = 1;
model.sim_type = 'NN';
model.mdl = {};

% create DFSM
model = createDFSM(data,model,iCase,state_names);

% for the second state run simulations using model from the first state
iCase = 2;
model = createDFSM(data,model,iCase,state_names);

end


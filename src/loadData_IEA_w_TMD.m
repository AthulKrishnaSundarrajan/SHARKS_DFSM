function data = loadData_IEA_w_TMD

close all;

root_path = which('INSTALL_DFSM'); % obtain full function path
data_path = fullfile(fileparts(root_path), 'data', filesep);
fulldata_path = fullfile(data_path,'DFSM_step'); % new data

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
req_states = {'PtfmPitch','TTDspFA','GenSpeed'};

req_controls = {'Wind1VelX','GenTq','BldPitch1'};
n_names = length(output_names);iTime = 1;
sim_plot = 0;

wind = cell(nLinCases,1);

reduce_samples = true;

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

        nt = 1000;
        t_ = linspace(t(1),t(end),nt);
        x = interp1(t,x,t_,'pchip');
        u = interp1(t,u,t_,'pchip');
        t = t_;



    end

    % assign
    data(iCase).time = t;
    data(iCase).states = x;
    data(iCase).inputs = u;
    data(iCase).state_names = ChanName(iStates);
    data(iCase).input_names = ChanName(iInputs);

    wind{iCase} = u(:,1);

   
end



% specify which states derivative is included
dindex = 1:length(req_states);

% approximate state derivatives
data = approximateStateDerivatives(data,dindex,1,1);

% add state names
LinearModels = [];
state_names = req_states;
for k = 1:length(dindex)
    state_names{end+1} = ['d',req_states{dindex(k)}];
end

% generate model using first simulation
t_case = 2;
model.sim_type = 'LTI';
model.mdl = {};

% create DFSM
model = createDFSM(data,model,t_case,state_names,req_controls);

% for the second state run simulations using model from the first state
v_case = 3;
model = createDFSM(data,model,v_case,state_names,req_controls);

% plot histogram
hf = figure;
hf.Color = 'w';
hold on;

histogram(wind{t_case})
histogram(wind{v_case})

legend('T','V')


tcase_range = [max(wind{t_case}),min(wind{t_case})];

hf = figure;
hf.Color = 'w';

hold on;

yline(tcase_range);
plot(data(v_case).time,wind{v_case});

end
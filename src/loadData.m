function data = loadData

% get current directory

current_folder = pwd;

% add DanZs openfast output processing toolbox to path
tool_path = [current_folder,'/matlab-toolbox'];
addpath(genpath(tool_path))

% path to folder where the openfast outputs are
Out_path = 'data';

% prefix and suffix of output files
prefix = 'IEA_w_TMD_';
suffix = '.outb';

% get all the files in the directory
Out_files = dir(fullfile(Out_path,[prefix,'*',suffix]));

% get length
nLinCases = length(Out_files);
% nLinCases = 1;

% numbering scheme bases on the number of files
if nLinCases <= 10
    numstring = '%01d';
else
    numstring = '%02d';
end

% reqired outputs
output_names = {'RtVAvgxh','GenTq','BldPitch1','PtfmPitch','GenSpeed','GenPwr'};
n_names = length(output_names);
sim_plot = 1;

% go through each case
for iCase = 1:nLinCases

    switch suffix
        case '.outb'
            [Channels, ChanName, ChanUnit, DescStr] = ReadFASTbinary(fullfile(Out_path,[prefix,num2str(iCase-1,'%01d'),suffix]));
        case '.mat'
            load([prefix,num2str(iCase-1,'%01d'),suffix]);
    end

    %go through each output name and plot trajectory
    if sim_plot
%         figure(iCase)
        for idx = 1:n_names
%             subplot(3,2,idx)
            ind = contains(ChanName,output_names{idx});
            %      find(ind)
%             plot(Channels(:,1),Channels(:,ind))
%             title([output_names{idx},'',ChanUnit{ind}])
%             xlim([100,800])
        end
    end
    iTime = 1;
    iStates = [23,9]; %,22,24,25,26,27]; % PtfmPitch, GenSpeed
%     iInputs = [2,268,6]; % RtVAvgxh, GenTq, BldPitch1
    % iInputs = [2,268]; % RtVAvgxh, GenTq, BldPitch1
    iInputs = [264,268]; % RtVAvgxh, GenTq, BldPitch1
    

    % extract
    t = Channels(:,iTime);
    x = Channels(:,iStates); x(:,1) = rad2deg(x(:,1));
    u = Channels(:,iInputs);

    I = 28; plot(t, Channels(:,I)); title(ChanName{I})

    % assign
    data(iCase).time = t;
    data(iCase).states = x;
    data(iCase).inputs = u;
    data(iCase).state_names = ChanName(iStates);
    data(iCase).input_names = ChanName(iInputs);

end

end
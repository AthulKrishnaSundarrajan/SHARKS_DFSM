function data = loadOFsim(fol_name,prefix,suffix,req_states,req_controls)

root_path = which('INSTALL_DFSM'); % obtain full function path
data_path = fullfile(fileparts(root_path), 'data', filesep);
fulldata_path = fullfile(data_path,fol_name); % new data

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

iTime = 1;


% go through each case
for iCase = 1:nLinCases

    switch suffix
        case '.outb'
            [Channels, ChanName, ChanUnit, DescStr] = ReadFASTbinary(fullfile(fulldata_path,[prefix,num2str(iCase-1,numstring),suffix]));
        case '.mat'
            load([prefix,num2str(iCase-1,'%01d'),suffix]);
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



end
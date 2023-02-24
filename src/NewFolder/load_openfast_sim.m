% function to load openfast simulation

function sim_detail = load_openfast_sim(sim_name,reqd_states,reqd_controls,tmin,tmax)

    % find type of binary file
    suffix = sim_name(end-4:end);
    
    % check the type of outputfile and load results
    if strcmpi(suffix,'.outb')

        % binary file (faster to read, and smaller size)

        [Channels, ChanName, ChanUnit, DescStr] = ReadFASTbinary(sim_name);

    elseif strcmpi(suffix,'.out')

        % text file(takes longer to read, and has a larger file size)

        [Channels, ChanName, ChanUnit, DescStr] = ReadFASTtext(sim_name);


    end
    
    % extract time
    time = Channels(:,1);
    time = time - min(time);

    if (tmin == 0) && isempty(tmax)

        t_ind = true(length(time),1);

    elseif ~isempty(tmax)

        t_ind = (time<=tmax) & (time>=tmin);

        time = time(t_ind);

    end
    
    % extract state index
    for i = 1:length(reqd_states)
        state_ind(i) = find(contains(ChanName,reqd_states{i}));
    end

    % extract control index
    for i = 1:length(reqd_controls)
        control_ind(i) = find(contains(ChanName,reqd_controls{i}));
    end

    % extract channels
    states = Channels(t_ind,state_ind);
    controls = Channels(t_ind,control_ind);

    % add to struct
    sim_detail.time = time;
    sim_detail.states = states;
    sim_detail.controls = controls;
    sim_detail.state_names = reqd_states;
    sim_detail.control_names = reqd_controls;


end
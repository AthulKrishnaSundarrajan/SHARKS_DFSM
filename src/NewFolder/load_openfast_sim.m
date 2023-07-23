% function to load openfast simulation

function sim_detail = load_openfast_sim(sim_name,reqd_states,reqd_controls,reqd_outputs,scale_outputs,tmin,tmax)

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
    ptfm_surge = false(length(reqd_states),1);
    for i = 1:length(reqd_states)
        state_ind(i) = find(contains(ChanName,reqd_states{i}));

        if strcmpi(reqd_states{i},'PtfmSurge')
            ptfm_surge(i) = true;
        end
    end
    
    if ~isempty(reqd_outputs)
        for i = 1:length(reqd_outputs)
            output_ind(i) = find(contains(ChanName,reqd_outputs{i}));
        end

        outputs = Channels(t_ind,output_ind);
        outputs = outputs;

        if scale_outputs

            scale_fac = 10.^(floor(log10(mean(outputs))));
            outputs = outputs./scale_fac;



        end

    else
        output_ind = [];
        outputs = [];

            
    end

    % extract control index
    gt_flag = false(length(reqd_controls),1);
    for i = 1:length(reqd_controls)
        control_ind(i) = find(contains(ChanName,reqd_controls{i}));

        if strcmpi(reqd_controls{i},'GenTq')
            gt_flag(i) = true;
        end
    end

    % extract channels
    states = Channels(t_ind,state_ind);
    controls = Channels(t_ind,control_ind);

    controls(:,gt_flag) = controls(:,gt_flag)/1000;
    states(:,ptfm_surge) = states(:,ptfm_surge)/1;

    % add to struct
    sim_detail.time = time;
    sim_detail.states = states;
    sim_detail.controls = controls;
    sim_detail.state_names = reqd_states;
    sim_detail.control_names = reqd_controls;
    sim_detail.outputs = outputs;
    sim_detail.output_names = reqd_outputs;


end
% function to filter signals
function sim_details = filter_signals(sim_details,filter_args)

    % extract details
    time = sim_details.time;
    states = sim_details.states;
    controls = sim_details.controls;
    outputs = sim_details.outputs;

    % extract arguments
    filt_states = filter_args.filt_states;
    filt_states_tf = filter_args.filt_states_tf;

    filt_controls = filter_args.filt_controls;
    filt_controls_tf = filter_args.filt_controls_tf;
    
    % filter states
    for ist = 1:length(filt_states)
        
        if filt_states(ist)

            t_f = filt_states_tf(ist);
            states(:,ist) = perform_filter(states(:,ist),time,t_f);

        end

    end
    
    % filter controls
    for ict = 1:length(filt_controls)
        
        if filt_controls(ict)

            t_f = filt_controls_tf(ict);
            controls(:,ict) = perform_filter(controls(:,ict),time,t_f);

        end

    end

    % filter outputs
    if ~isempty(outputs)

         filt_outputs = filter_args.filt_outputs;
         filt_outputs_tf = filter_args.filt_outputs_tf;


        for iot = 1:length(filt_outputs)

            if filt_outputs(iot)

                t_f = filt_outputs_tf(iot);
                outputs(:,iot) = perform_filter(outputs(:,iot),time,t_f);

            end
        end


    end

    % add

    sim_details.states = states;
    sim_details.controls = controls;
    sim_details.outputs = outputs;


end

function signal = perform_filter(signal,time,t_f)

    if t_f == 0
    
    else

    dt = time(2) - time(1);

    %
    nb = floor(t_f/dt);
    
    %
    b = ones(nb,1)/nb;
    
    
    %
    signal = filtfilt(b,1,signal);
    end

end
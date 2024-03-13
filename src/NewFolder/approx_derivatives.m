% function to approximate derivatives
function sim_details = approx_derivatives(sim_details,add_dx2)

    % extract time and state information
    time = sim_details.time;
    states = sim_details.states;

    % construct polynomial approximation
    states_pp = spline(time,states');
    
    % construct derivative
    dx_pp = fnder(states_pp,1);

    dx2_pp = fnder(states_pp,2);
    
    % evaluate state derivative
    state_derivatives = ppval(dx_pp,time);

    state_derivatives2 = ppval(dx2_pp,time);

    if size(states,2) == 1
        state_derivatives = state_derivatives';
        state_derivatives2 = state_derivatives2';
    end
    
    % assign
    if add_dx2
        states = [states,state_derivatives'];
        state_derivatives = [state_derivatives;state_derivatives2];
        sim_details.states = states;

        state_names = sim_details.state_names;
        dx_names = cell(size(state_names));
        for idx  = 1:length(state_names)
            dx_names{idx} = ['d',state_names{idx}];
        end
% 
         sim_details.state_names = horzcat(state_names,dx_names);
    end

    sim_details.state_derivatives = state_derivatives';

end
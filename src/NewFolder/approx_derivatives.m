% function to approximate derivatives
function sim_details = approx_derivatives(sim_details)

    % extract time and state information
    time = sim_details.time;
    states = sim_details.states;

    % construct polynomial approximation
    states_pp = spline(time,states');
    
    % construct derivative
    dx_pp = fnder(states_pp,1);
    
    % evaluate state derivative
    state_derivatives = ppval(dx_pp,time);
    
    % assign
    sim_details.state_derivatives = state_derivatives';

end
function sim_details = cell2struct_dfsm(sim_details_cell,nsim)

% get length
sim_details(nsim) = struct();

% loop through
for i = 1:nsim
    
    % initialize
    sim_details(i).time = [];sim_details(i).states = [];sim_details(i).controls = [];
    sim_details(i).state_names = [];sim_details(i).control_names = [];sim_details(i).state_derivatives = [];
    sim_details(i).nstates = [];sim_details(i).ncontrols = [];sim_details(i).ninputs = [];sim_details(i).noutputs = [];
    sim_details(i).dx_act = [];sim_details(i).fun_name = [];

    % assign
    sim_details(i) = sim_details_cell{i};

end

end
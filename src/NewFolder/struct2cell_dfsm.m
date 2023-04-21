function [inputs_cell,state_dx_cell,outputs_cell] = struct2cell_dfsm(sim_details)

    % get the number of samples
    nsamples = length(sim_details);
    
    inputs_cell = cell(nsamples,1);
    state_dx_cell = cell(nsamples,1);
    outputs_cell = cell(nsamples,1);
    
    % loop through and extract cell
    for isample = 1:nsamples
    
        controls = sim_details(isample).controls;
        states = sim_details(isample).states;
        state_derivatives = sim_details(isample).state_derivatives;
        outputs = sim_details(isample).outputs;
    
        % combine
        inputs = [controls,states];
        
        % store
        inputs_cell{isample} = inputs;
        state_dx_cell{isample} = state_derivatives;
        outputs_cell{isample} = outputs;
    
    end
end
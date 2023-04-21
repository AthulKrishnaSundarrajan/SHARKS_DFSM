function [input_sampled,state_dx_sampled,output_sampled,inputs,state_dx,outputs] = sample_data(train_samples,sampling_type,n_samples)

    % function to sample the inputs and outputs
    % two sampling stratergies available:
    %   1. equi-distant sampling - 'ED'
    %   2. k-means clustering sampling - 'KM'
    
    % extract data from structure
    if isstruct(train_samples)
        [inputs,state_dx,outputs] = struct2cell_dfsm(train_samples);

        % concatenate
        inputs = vertcat(inputs{:});
        state_dx = vertcat(state_dx{:});
        outputs = vertcat(outputs{:});

    elseif iscell(train_samples)

        inputs = train_samples{1};
        state_dx = train_samples{2};
        outputs = train_samples{3};


    end

    if isnan(n_samples)
        input_sampled = inputs;
        state_dx_sampled = state_dx;
        output_sampled = outputs;

        return
    end


    % sample
    switch sampling_type
 
    
        case 'KM'
            % perform k means clustering
            [input_sampled,state_dx_sampled,output_sampled] = perform_KM(inputs,state_dx,outputs,n_samples);
    
        case 'random'
    
            % get number of inputs
            nt = size(inputs,1);
            
            % generate random
            index = randsample(nt,n_samples);
    
            % extract
            input_sampled = inputs(index,:);
            state_dx_sampled = state_dx(index,:);
            output_sampled = outputs(index,:);

    end


end
function[input_sampled,state_dx_sampled,output_sampled] = perform_KM(inputs,state_dx,outputs,n_samples)

    % get clusters
    [idx,input_sampled] = kmeans(inputs,n_samples,'MaxIter',1000);

    state_dx_sampled = zeros(n_samples,size(state_dx,2));
    output_sampled = zeros(n_samples,size(outputs,2));

    for icluster = 1:n_samples
        
        % find the index corresponding to given cluster
        ind_cluster = (idx == icluster);
        
        % find the state derivative values for those clusters
        state_dx_cluster = state_dx(ind_cluster,:);

        if ~isempty(outputs)
         output_cluster = outputs(ind_cluster,:);

         output_sampled(icluster,:) = mean(output_cluster,1);
        end
        
        % the state derivative value at the cluster centroid is the
        % mean
        state_dx_sampled(icluster,:) = mean(state_dx_cluster,1);

    end

end
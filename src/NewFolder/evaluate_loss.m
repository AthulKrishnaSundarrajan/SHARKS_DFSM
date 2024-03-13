function V = evaluate_loss(x,ns,nc,inputs,dx_act)

    % get linear model
    [A,B,C,D] = LTI_function(x,[],ns,nc);
    
    % evaluate state derivatives
    dx_predicted = inputs*[B,A]';

    % number of data samples
    N = length(inputs);

    % calculate error
    error = dx_act - dx_predicted;
  
    % calculate loss
    V = 1/N*(trace(error'*error));


end

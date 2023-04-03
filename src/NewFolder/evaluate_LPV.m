function dx = evaluate_LPV(lin,wmax,wmin,nstates,inputs,noutputs)

% function to evaluate the LPV model for a given set of inputs
% works for both direct evaluation and through ode

% get the size of inputs
nt = size(inputs,1);

% initialize derivative
dx = zeros(nt,noutputs);

% go through the inputs and evaluate lpv model
for i = 1:nt

    % extract wind
    w = inputs(i,1);
    
    % check value to see if the wind speed is between limits
    if w > wmax
        w = wmax;
    elseif w<wmin
        w = wmin;
    end
    
    % evaluate derivative function
    dx(i,:) = inputs(i,:)*lin(w);
end

end
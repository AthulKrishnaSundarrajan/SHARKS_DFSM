function V = grey_box_objective(x,ns,nc,inputs,dx_act,time,sim_flag)


% Calculate loss
[V,A,B,C,D] = evaluate_loss(x,ns,nc,inputs,dx_act);

n_exp = length(inputs)/(length(time));

% if sim_flag is set to true, then a part of the data is simulated using
% the surrogate model

if sim_flag
    V_sim = evaluate_loss_sim(A,B,inputs,time,n_exp,nc,ns);
    
    V = V_sim;
end
end



function [V,A,B,C,D] = evaluate_loss(x,ns,nc,inputs,dx_act)

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

function V_sim = evaluate_loss_sim(A,B,inputs,time,n_exp,nc,ns)

% get the number of points
nt = length(time);

% value to scale states by
states_max = max(abs(inputs(:,nc+1:end)));
states_max = states_max(1:ns/2);

% reshape inputs from [(nt*n_exp)x(nc+ns)] to [(n_exp)x(nt)x(nc+ns)]
inputs_re = reshape(inputs,[n_exp,nt,nc+ns]);

% initalize storage array
loss_exp = zeros(n_exp,1);

% loop through and simulate a part of the data to evaluate the loss between
% the actual states and predicted states
for i = 1:n_exp
    
    % get the inputs for the current experiment
    inputs_ = squeeze(inputs_re(i,:,:));
    
    % extract controls and states
    controls = inputs_(:,1:nc);
    states = inputs_(:,nc+1:end);

    % get t_max
    t_max = max(time);

    % simulate the last 100 secs
    t_start_test = t_max - 500;
    
    % get the index for the data points between t_start_test and t_max
    t_ind = (time>=t_start_test) & (time <= t_max);
    
    % extract states and controls corresponding to t_ind
    states_test = states(t_ind,:);
    controls_test = controls(t_ind,:);
    time_test = time(t_ind);
    
    % output matrices
    C = [eye(ns/2),zeros(ns/2)];
    D = zeros(ns/2,nc);
    
    % construct state space
    sys = ss(A,B,C,D);
    
    % initial states
    X0 = states_test(1,:);
    
    % simulate the system using lsim
    states_pred = lsim(sys,controls_test,time_test,X0,'zoh');

    % scale
    states_act = states_test(:,1:ns/2);
    states_test = states_act./states_max;
    states_pred = states_pred./states_max;

    % evaluate the error
    error = states_test - states_pred;
    
    % calculate loss
    loss_exp(i) = 1/(length(time_test))*(trace(error'*error));

end

% mean loss
V_sim = mean(loss_exp);
end


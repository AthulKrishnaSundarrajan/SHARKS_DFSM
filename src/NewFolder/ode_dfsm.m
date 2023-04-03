function dx = ode_dfsm(t,x,u_fun,dfsm)
    
    % controls
    u = u_fun(t);
    u = u(:);

    % combine
    inputs = [u',x'];
    
    % evaluate dfsm
    dx = evaluate_dfsm(inputs,dfsm,'deriv');

end
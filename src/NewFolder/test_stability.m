function [c,ceq] =  test_stability(x,ns,nc)

    % get the linear model
    [A,~,~,~] = LTI_function(x,[],ns,nc);
    
    % evaluate the eigen values
    eigA = eig(A);
    
    % real values
    c = real(eigA);
    
    % no equality constraints
    ceq = [];

end
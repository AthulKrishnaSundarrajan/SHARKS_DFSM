function f = construct_function(AB,nonlin,error_ind,fun_type,nx,nu)

% initialize
f = cell(nx,1);

% calculate the number of inputs
ninputs = nx+nu;

% loop through and construct anonymous functions
for ix = 1:nx

    % extract
    nx_fun = nonlin{ix};
    
    % construct function
    if error_ind(ix)

        switch fun_type

            case 'GPR'

                switch class(AB)

                    case 'function_handle'

                            f{ix} = @(t,param,UYP) 

                    case 'double'
                        if ~isempty(AB)
        
                            f{ix} = @(t,param,UYP) UYP(:,1:ninputs)*AB(:,ix) + predict(nx_fun,UYP(:,1:ninputs));
        
                        else
        
                            f{ix} = @(t,param,UYP) predict(nx_fun,UYP(:,1:ninputs));
        
                        end

                end

            case 'RBF'

                if ~isempty(AB)

                    f{ix} = @(t,param,UYP) UYP(:,1:ninputs)*AB(:,ix) + (nx_fun((UYP(:,1:ninputs))')');

                else

                    f{ix} = @(t,param,UYP)(nx_fun((UYP(:,1:ninputs))')');
                end
        end

    else 
        f{ix} = @(t,param,UYP) UYP(:,1:ninputs)*AB(:,ix);
    end


end

end
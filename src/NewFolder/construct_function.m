function f = construct_function(AB,nonlin,error_ind,fun_type,type,nx,nu,ny,con)

% initialize

if strcmpi(type,'deriv')
    f = cell(nx,1);

    con_val = zeros(nx,1);
elseif strcmpi(type ,'op')
    f = cell(ny,1);

    con_val = con;
end

nf = length(f);

% calculate the number of inputs
ninputs = nx+nu;

% loop through and construct anonymous functions
for ind_f = 1:nf

    % extract
    nx_fun = nonlin{ind_f};
    
    % construct function
    if error_ind(ind_f)

        switch fun_type

            case 'GPR'

                switch class(AB)

                    case 'function_handle'

                            %f{ix} = @(t,param,UYP) 

                    case 'double'
                        if ~isempty(AB)
        
                            f{ind_f} = @(t,param,UYP) UYP(:,1:ninputs)*AB(:,ind_f) + predict(nx_fun,UYP(:,1:ninputs)) - con_val(ind_f);
        
                        else
        
                            f{ind_f} = @(t,param,UYP) predict(nx_fun,UYP(:,1:ninputs))- con_val(ind_f);
        
                        end

                end

            case 'RBF'

                if ~isempty(AB)

                    f{ind_f} = @(t,param,UYP) UYP(:,1:ninputs)*AB(:,ind_f) + (nx_fun((UYP(:,1:ninputs))')')- con_val(ind_f);

                else

                    f{ind_f} = @(t,param,UYP)(nx_fun((UYP(:,1:ninputs))')')- con_val(ind_f);
                end
        end

    else 
        f{ind_f} = @(t,param,UYP) UYP(:,1:ninputs)*AB(:,ind_f)- con_val(ind_f);
    end


end

end
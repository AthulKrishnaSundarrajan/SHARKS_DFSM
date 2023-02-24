clc;clear; close all;







function dx = og_deriv_function(t,x,u_fun,fun_name)

    % control function
    u = u_fun(t);
    

    switch fun_name

        case 'vanderpol'
            
            % extract
            x1 = x(1); x2 = x(2);

            % calculate derivative
            dx = [x2;-x1 + x2 - x1^2*x2 + u];

        case 'math1'
            
            % extract
            x1 = x(1); x2 = x(2);

            % calculate state derivative
            a = 0.7137;b = 0.3929;
            dx = [-a*x1 + b*x2^2;b*x1 - 2*a^2*x2 - x1^2 + u];

        case 'math2'
            
            % extract
            x1 = x(1);x2 = x(2);

            % calculate state derivative
            a = -0.5204;b = 0.4926;
            dx = [-0.84*x1 - a*x2 -b*x1*x2; 0.54*x1 + a*x2 + b*x1*x2 + u];

    end


end
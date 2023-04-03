% some example nonlinear functions

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
            dx = [-a*x1+b*x2^2;
                b*x1-2*a^2*x2-x1^2+u];

        case 'math2'
            
            % extract
            x1 = x(1);x2 = x(2);

            % calculate state derivative
            a = -0.5204;b = 0.4926;
            dx = [-0.84*x1-a*x2-b*x1*x2;
                0.54*x1+a*x2+b*x1*x2+u];

        case 'two-link-robot'

            % extract
            x1 = x(1);x2 = x(2);x3 = x(3);
            
            % controls
            u1 = u(1);u2 = u(2);
            
            % calculate state derivative
            dx = [(sin(x3)*(9/4*cos(x3)*x1^2+2*x2^2)+4/3*(u1-u2)-3/2*cos(x3)*u2)/(31/36+9/4*sin(x3)*sin(x3));
                -(sin(x3)*(7/2*x1^2+9/4*cos(x3)*x2^2)-7/3*u2+3/2*cos(x3)*(u1-u2))/(31/36+9/4*sin(x3)*sin(x3));
                x2-x1;
                x1];

        case 'transfer-min-fuel'

            % extract
            x1 = x(1);x2 = x(2);x3 = x(3);x4 = x(4);

            % controls
            u1 = u(1);u2 = u(2);

            % calculate state derivatives
            dx = [x3;
                  x4/x1;
                  x4^2/x1-1/x1^2+u1;
                  -x4*x3/x1+u2];

    end


end
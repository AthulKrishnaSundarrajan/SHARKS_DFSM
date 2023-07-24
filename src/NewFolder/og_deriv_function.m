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

        case 'two-link-robot2'

            J1 = 0.12;
            J2 = 0.15;
            g = 9.81;
            l1 = 0.6;
            l2 = 0.8;
            m1 = 3;
            m2 = 2.5;

            % extract state and control variables
            tau1 = u(1); tau2 = u(2);
            T1 = x(1); DT1 = x(2); T2 = x(3); DT2 = x(4);
        
            % extract system parameters
%             J1 = p.J1; J2 = p.J2; g = p.g; l1 = p.l1; l2 = p.l2;
%             m1 = p.m1; m2 = p.m2;
        
        
            dx2 = (4*tau1*(m2*l2^2 + 4*J2))/(16*J1*J2 + 4*l1^2*l2^2*m2^2 + 4*J2*l1^2*m1 + 4*J1*l2^2*m2 + 16*J2*l1^2*m2 + l1^2*l2^2*m1*m2- 4*l1^2*l2^2*m2^2*cos(T2)^2) - (4*tau2*(m2*l2^2 + 2*l1*m2*cos(T2)*l2 + 4*J2))/(16*J1*J2 + 4*l1^2*l2^2*m2^2 + 4*J2*l1^2*m1+ 4*J1*l2^2*m2 + 16*J2*l1^2*m2 + l1^2*l2^2*m1*m2 - 4*l1^2*l2^2*m2^2*cos(T2)^2) + (4*((l1*l2*m2*sin(T2)*DT1^2)/2+ (g*l2*m2*cos(T1 + T2))/2)*(m2*l2^2 + 2*l1*m2*cos(T2)*l2 + 4*J2))/(16*J1*J2 + 4*l1^2*l2^2*m2^2 + 4*J2*l1^2*m1 + 4*J1*l2^2*m2 + 16*J2*l1^2*m2+ l1^2*l2^2*m1*m2 - 4*l1^2*l2^2*m2^2*cos(T2)^2) - (4*(m2*l2^2 + 4*J2)*(g*m2*((l2*cos(T1 + T2))/2 + l1*cos(T1)) + (g*l1*m1*cos(T1))/2 - (DT2*l1*l2*m2*sin(T2)*(2*DT1 + DT2))/2))/(16*J1*J2 + 4*l1^2*l2^2*m2^2 + 4*J2*l1^2*m1 + 4*J1*l2^2*m2 + 16*J2*l1^2*m2 + l1^2*l2^2*m1*m2 - 4*l1^2*l2^2*m2^2*cos(T2)^2);
            
            dx4 = (4*(g*m2*((l2*cos(T1 + T2))/2 + l1*cos(T1)) + (g*l1*m1*cos(T1))/2 ...
                - (DT2*l1*l2*m2*sin(T2)*(2*DT1 + DT2))/2)*(m2*l2^2 + 2*l1*m2*cos(T2)*l2 + 4*J2))/(16*J1*J2 + 4*l1^2*l2^2*m2^2 + 4*J2*l1^2*m1 + 4*J1*l2^2*m2 + 16*J2*l1^2*m2 + l1^2*l2^2*m1*m2 - 4*l1^2*l2^2*m2^2*cos(T2)^2)...
                + (4*tau2*(4*J1 + 4*J2 + l1^2*m1 + 4*l1^2*m2 + l2^2*m2 + 4*l1*l2*m2*cos(T2)))/(16*J1*J2 + 4*l1^2*l2^2*m2^2 + 4*J2*l1^2*m1 + 4*J1*l2^2*m2 + 16*J2*l1^2*m2 + l1^2*l2^2*m1*m2 - 4*l1^2*l2^2*m2^2*cos(T2)^2)...
                - (4*((l1*l2*m2*sin(T2)*DT1^2)/2 + (g*l2*m2*cos(T1 + T2))/2)*(4*J1 + 4*J2 + l1^2*m1 + 4*l1^2*m2 + l2^2*m2 + 4*l1*l2*m2*cos(T2)))/(16*J1*J2 + 4*l1^2*l2^2*m2^2 + 4*J2*l1^2*m1 + 4*J1*l2^2*m2 + 16*J2*l1^2*m2 + l1^2*l2^2*m1*m2 - 4*l1^2*l2^2*m2^2*cos(T2)^2) -...
                (4*tau1*(m2*l2^2 + 2*l1*m2*cos(T2)*l2 + 4*J2))/(16*J1*J2 + 4*l1^2*l2^2*m2^2 + 4*J2*l1^2*m1 + 4*J1*l2^2*m2 + 16*J2*l1^2*m2 + l1^2*l2^2*m1*m2 - 4*l1^2*l2^2*m2^2*cos(T2)^2);
        
        
            dx = [DT1;
                dx2;
                DT2;
                dx4];


    end


end
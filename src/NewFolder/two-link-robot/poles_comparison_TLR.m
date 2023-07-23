clc; clear; close all;

%%

% Script to setup and run the frequency domain validation study from Sec 3
% The example discussed is the two-link robotic (TLR) system from Ex 2.10 in http://doi.org/10.1007/978-3-540-78486-9_4

%% Define the nonlinear system

% system parameters
p.J1 = 0.12;
p.J2 = 0.15;
p.g = 9.81;
p.l1 = 0.6;
p.l2 = 0.8;
p.m1 = 3;
p.m2 = 2.5;

% test
u_fun = @(t)[1,1]; %[26.1239,9.4757];
x = [1,1,1,1]';% [deg2rad(45);0;deg2rad(-30);0]; 
t = 1;

dx = TLR_nonlin(t,x,u_fun,p);

% linearized matrices
A_taylor =  [0,1.0000,0,0;
    17.8320,0,-3.3656,0
    0,0,0,1.0000;
  -30.0624,0,16.4148,0];

B_taylor = [0,0;1.2514,-2.4337;0,0;-2.4337,6.5512];

C = [1 0 0 0 ; 0 0 1 0];
D = zeros(2);

return
%--------------------------------------------------------------------------
function dx = TLR_nonlin(t,x,u_fun,p)

    %% nonlinear derivative function for the TLR example
    
    % evaluate control inputs
    u = u_fun(t);

    % extract state and control variables
    tau1 = u(1); tau2 = u(2);
    T1 = x(1); DT1 = x(2); T2 = x(3); DT2 = x(4);

    % extract system parameters
    J1 = p.J1; J2 = p.J2; g = p.g; l1 = p.l1; l2 = p.l2;
    m1 = p.m1; m2 = p.m2;

    val = [x;u';J1;J2;g;l1;l2;m1;m2];

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

function dx_taylor = test_taylor(A,B,test_simulations,x0,u0)

    ntest = length(test_simulations);

    dx_taylor = cell(ntest,2);

    for i = 1:ntest

        states = test_simulations(i).states;
        controls = test_simulations(i).controls;
        state_derivatives = test_simulations(i).state_derivatives;

        dx = A*(states-x0' ) + B*(controls-u0');

        dx_taylor{i,1} = state_derivatives;
        dx_taylor{i,2} = dx;
        

    end
    
end

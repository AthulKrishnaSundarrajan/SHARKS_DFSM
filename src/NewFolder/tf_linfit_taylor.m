clc; clear; close all;

% set seed
rng(4357)

% define parameters
t0 = 0; t_f = 0.1;
fname = mfilename('fullpath');


fname = which(fname);

nt = 90; nsamples = 100;
fun_name = 'two-link-robot2';
split = [1,0];

% run simulation and get results
nfac = 5;
fac_ = [0,-1];
nfac = length(fac_);
fac = 10.^fac_;

% run simulation and get results

%fac = 1e-6; %logspace(-2,0,nfac);
x0 = [deg2rad(45),0,deg2rad(-30),0];
u0 = [26.1239,9.4757];

dfsm_options.ltype = 'LTI';
dfsm_options.ntype = '';
dfsm_options.nsamples = nan;
dfsm_options.sampling_type = 'KM';
dfsm_options.train_test_split = split;
dfsm_options.scale_flag = ~true;


%sim_details = run_simulation(t0,t_f,nt,nsamples,fun_name,fac,x0,u0);


saveflag = false;
fol_name = 'plots_linear_validation';
x_lim = [0,t_f];


%-------------------------------------------------
% define symbolic variables
% syms x1 x2 x3 x4 u1 u2 X1 X2 X3 X4 U1 U2 real
% 
% % derivative function
% x1_d = ((sin(x3).*(9/4*cos(x3).*x1.^2+2*x2.^2) + 4/3*(u1-u2) - 3/2*cos(x3).*u2 )./ (31/36 + 9/4*sin(x3).^2) );
% x2_d = (-( sin(x3).*(9/4*cos(x3).*x2.^2+7/2*x1.^2) - 7/3*u2 + 3/2*cos(x3).*(u1-u2) )./ (31/36 + 9/4*sin(x3).^2) );
% x3_d = ( x2-x1 ); 
% x4_d = ( x1 );
% 
% % outputs
% dx = [x1_d;x2_d;x3_d;x4_d];
% 
% % inputs
% x = [x1;x2;x3;x4];
% u = [u1;u2];
% 
% %x0 = [0;0;0.5;0];
% %u0 = [0;0];
% 
% I0 = [u0;x0];
% 
% 
% % inputs
% I = [u;x];
% I_ = [U1;U2;X1;X2;X3;X4];
% 
% % evaluate jacobian
% L_taylor = jacobian(dx,I);
% 
% % substitute value
% L_taylor = double(subs(L_taylor,I,I0));
% 
% % extract
% A_taylor = L_taylor(:,3:end);
% B_taylor = L_taylor(:,1:2);
% 

A_taylor =  [0,1.0000,0,0;
    17.8320,0,-3.3656,0
    0,0,0,1.0000;
  -30.0624,0,16.4148,0];

% A_taylor =  [0,1.0000,0,0;
%     17.8320,0,-3.0024,0
%     0,0,0,1.0000;
%   -30.0624,0,10.456,0];


B_taylor = [0,0;1.2514,-2.4337;0,0;-2.4337,6.5512];

C_ = [1 0 0 0 ; 0 0 1 0];
D = zeros(2);


sys_taylor = ss(A_taylor,B_taylor,C_,D);
tf_taylor = tf(sys_taylor);

linfit_cell = cell(nfac,1);



for i = 1:nfac
            sim_details = run_simulation(t0,t_f,nt,nsamples,fun_name,fac(i),x0,u0);
            
            train_simulations_ind = 1:floor(nsamples*0.8);
            train_simulations = sim_details(train_simulations_ind);
            
            test_simulations_ind = floor(nsamples*0.8+1):nsamples;
            test_simulations = sim_details(test_simulations_ind);
            
            dfsm = DFSM(train_simulations,dfsm_options);
            
            L_linfit = (dfsm.deriv.AB)';
            A_linfit = L_linfit(:,3:end);
            B_linfit = L_linfit(:,1:2);
            
            sys_linfit = ss(A_linfit,B_linfit,C_,D);
            
            tf_linfit = tf(sys_linfit);

            linfit_cell{i} = tf_linfit;      
        
end


%-----------------------------------------------------------

ind_pole = 3;
for ix = 1

    for iy = 1

     
        [z_taylor,p_taylor,~] = zpkdata(tf_taylor(ix,iy));
        p_taylor = p_taylor{1};
        
      
       
        hf = figure;
        hf.Color = 'w'; hold on;
        commonFigureProperties
        markersize = 30;
        
        xlabel('Real');ylabel('Imaginary')
        
        plot(real(p_taylor),imag(p_taylor),'k.','markersize',markersize)
        legend_entries = {'Taylor'};
        
        for i = 1:nfac

            linfit_ = linfit_cell{i};
            [z_linfit,p_linfit,~] = zpkdata(linfit_(ix,iy));
        
            legend_entries{end+1} = ['$\delta = 10^{',num2str(fac_(i)),'}$'];
            
            p_linfit_xy = p_linfit{1};
            
        
            plot(real(p_linfit_xy),imag(p_linfit_xy),'.','markersize',markersize)
        
        end
        
        %ylim([-0.007,0.007])
        taskflag = 'axes'; commonFigureTasks;
        %xlim([-0.04,0.04]);ylim([-0.01,0.01])

    end
end

legend(legend_entries)

fontlegend = 15; nCol = 3;lcn = 'northoutside';
taskflag = 'legend';commonFigureTasks;


if saveflag
    savename = 'poles_comparison';
    pathpdf = mfoldername(mfilename('fullpath'),fol_name);
    filename = fullfile(pathpdf,savename);
    str = strcat("export_fig '",filename,"' -pdf");
    eval(str)
end

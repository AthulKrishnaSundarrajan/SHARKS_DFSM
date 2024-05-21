clc; clear; close all;

rng(34534)

% file name

root_path = fileparts(which('INSTALL_DFSM'));
data_path = fullfile(root_path,'data');
fol_name_cell = 'transition';
fol_name2 = 'transition_test';


sim_path = fullfile(data_path,fol_name_cell);
sim_path2 = fullfile(data_path,fol_name2);

% file names
prefix = 'IEA_w_TMD_';
suffix = '.outb';

% sim_files 
sim_dir  = dir(fullfile(sim_path,[prefix,'*',suffix]));
sim_dir2 = dir(fullfile(sim_path2,[prefix,'*',suffix]));

nsim = length(sim_dir);

if nsim <= 10
    numstring = '%01d';
else
    numstring = '%02d';
end

sim_files = cell(nsim,1);

% store the name of the files
for i = 1:nsim
    sim_files{i} = fullfile(sim_path,[prefix,num2str(i-1,numstring),suffix]);
end

nsim = length(sim_dir2);

if nsim <= 10
    numstring = '%01d';
else
    numstring = '%02d';
end

sim_files2 = cell(nsim,1);

% store the name of the files
for i = 1:nsim
    sim_files2{i} = fullfile(sim_path2,[prefix,num2str(i-1,numstring),suffix]);
end



% required states
reqd_states = {'PtfmPitch','TTDspFA','GenSpeed'};% ,'GenSpeed','YawBrTAxp'

% filtering arguments for states
filter_args.filt_states = [true,true,~false];
filter_args.filt_states_tf = 0.1*[1,1,1]; % 2,2 works


% required controls
reqd_controls = {'RtVAvgxh','GenTq','BldPitch1','Wave1Elev'};
reqd_outputs = {'TwrBsMxt','TwrBsMyt','TwrBsMzt','TwrBsFxt','NcIMURAys','YawBrTAxp'};%{'TwrBsFxt','TwrBsMyt','GenPwr','YawBrTAxp'};
filter_args.filt_outputs = false*ones(length(reqd_outputs),1);
filter_args.filt_outputs_tf = zeros(length(reqd_outputs),1);


% filtering arguments for controls
filter_args.filt_controls = [~true,false,false];
filter_args.filt_controls_tf = [0,0,0];

% filter flag
filter_flag = ~false;

%reqd_outputs = {};

% time
tmin = 00;
tmax = [];
add_dx2 = true;

% extract
sim_details = simulation_details(sim_files,reqd_states,reqd_controls,reqd_outputs,false,filter_flag,filter_args,add_dx2,tmin,tmax);
%sim_details = reorganize_data(sim_details);

sim_details2 = simulation_details(sim_files2,reqd_states,reqd_controls,reqd_outputs,false,filter_flag,filter_args,add_dx2,tmin,tmax);
%sim_details2 = reorganize_data(sim_details2);

train_index =  {1:5,8:12,15:19,22:26,29:33};%{1:5,10:15,20:25}; %%{1:5,10:15,20:25}
test_index = {6,13,20,27,34};%{}

nw = length(train_index);
x_cell = cell(nw,1);
A_cell = cell(nw,1);
B_cell = cell(nw,1);
C_cell = cell(nw,1);
D_cell = cell(nw,1);
model_construct_time = zeros(nw,1);

for iw = 1:nw

    nc = length(reqd_controls);
    ind_ = train_index{iw};
    [~,~,~,inputs,state_dx,outputs] = sample_data(sim_details(ind_),'KM',5);
    AB = linsolve(inputs,state_dx)';

    CD = linsolve(inputs,outputs);
    CD = CD';
    C = CD(:,nc+1:end);
    D = CD(:,1:nc);
    C_cell{iw} = C;
    D_cell{iw} = D;
    
    % 
    % % initial parameter
    ns = 2*length(reqd_states);
    par0 = AB(ns/2+1:end,:);
    
    x0 = par0(:);
    nx = length(x0);
    
    %
    time = sim_details(1).time; dt = time(2) - time(1);
    controls = inputs(:,1:nc);
    
    nt = length(time);
    
    %% testing out GA 
    
    global_options_list = {'GA','multi-start','global-search'};
    
    global_options = 'GA';
    
    %sim_flag = true;
    options_fmincon = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxIterations',500,'MaxFunctionEvaluations',10000,'UseParallel',true,"EnableFeasibilityMode",true,'ConstraintTolerance',1e-3);
    tic
    
    switch global_options
    
        case 'GA'
            
            if iw == 1
            
                % ga options
                gaopt = optimoptions('ga','UseParallel',true,'Display','iter','MaxGenerations',nx*2);%,'HybridFcn',{@fmincon,options_hybrid},'PopulationSize',500);
                
                % find good starting point
                sim_flag = false;
                x_global = ga(@(x) grey_box_objective(x,ns,nc,inputs,state_dx,time,sim_flag),nx,[],[],[],[],[],[],@(x)test_stability(x,ns,nc),[],gaopt) ;
                x = fmincon(@(x) grey_box_objective(x,ns,nc,inputs,state_dx,time,sim_flag),x_global,[],[],[],[],[],[],@(x)test_stability(x,ns,nc),options_fmincon);
            
            else

                x = fmincon(@(x) grey_box_objective(x,ns,nc,inputs,state_dx,time,sim_flag),x_cell{iw-1},[],[],[],[],[],[],@(x)test_stability(x,ns,nc),options_fmincon);
      
            
            end
            x_cell{iw} = x;
        case 'multi-start'
    
            problem = createOptimProblem('fmincon','objective',@(x) grey_box_objective(x,ns,nc,inputs,state_dx,time,sim_flag),'x0',par0,'nonlcon',@(x)test_stability(x,ns,nc),'options',options_fmincon);
    
            method = MultiStart;
    
            [x,f] = run(method,problem,50);
    
          
    
        case 'global-search'
    
            problem = createOptimProblem('fmincon','objective',@(x) grey_box_objective(x,ns,nc,inputs,state_dx,time,sim_flag),'x0',par0,'nonlcon',@(x)test_stability(x,ns,nc),'options',options_fmincon);
    
            method = GlobalSearch;
    
            [x,f] = run(method,problem);
    
    
    end
    
    
    model_construct_time(iw) = toc;
    
    [A,B,C,D] = LTI_function(x,[],ns,nc);
    A_cell{iw} = A;
    B_cell{iw} = B;
    
    
    % for ind_test = test_index{iw}
    % 
    % 
    %     % time and control inputs
    %     time = sim_details(ind_test).time;
    %     controls = sim_details(ind_test).controls;
    % 
    %     dfsm.ltype = 'LTI';
    %     dfsm.ntype = '';
    % 
    %     dfsm.deriv.AB = [B,A]';
    %     dfsm.deriv.nonlin = [];
    %     dfsm.nderiv = ns;
    %     dfsm.deriv.error_ind = [];
    %     dfsm.scaler_input = ones(1,nc+ns);
    %     sim_details(ind_test).outputs = [];
    % 
    %     % test the dfsm
    %     [dfsm_mf,X_cell,dx_cell,Y_cell] = test_dfsm(dfsm,sim_details(ind_test),ind_test,true,true);
    % end

end

hf = figure;
hf.Color = 'w';

hold on;
xlabel('Real');ylabel('Imag')
W = [9,10,11,12,13];%[16,17,18] ;%
nw = length(W);

n = 1;

for i = 1:nw

    A_ = A_cell{i};

    eig_A = eig(A_);
    
    %eig_A = eig_A(n);
    plot(real(eig_A),imag(eig_A),'.','MarkerSize',20,'color','r')
    text(real(eig_A)+0.0,imag(eig_A),num2str(W(i)))



end

n_out = length(reqd_outputs);

A = zeros(nw,ns,ns);
B = zeros(nw,ns,nc);
C = zeros(nw,n_out,ns);
D = zeros(nw,n_out,nc);

for i = 1:nw
    A(i,:,:) = A_cell{i};
    B(i,:,:) = B_cell{i};

    C(i,:,:) = C_cell{i};
    D(i,:,:) = D_cell{i};

end

A_pp = interp1(W,A,'nearest','pp');
A_fun = @(w) ppval(A_pp,w);

B_pp = interp1(W,B,'nearest','pp');
B_fun = @(w) ppval(B_pp,w);

C_pp = interp1(W,C,'nearest','pp');
C_fun = @(w) ppval(C_pp,w);

D_pp = interp1(W,D,'nearest','pp');
D_fun = @(w) ppval(D_pp,w);


test_inds = [1,4,7,10];%[3,8];

test_time = [0];
for i_ind = test_inds
    [~,~,~,inputs,state_dx,outputs] = sample_data(sim_details2(i_ind),'KM',5);
    time = sim_details2(1).time;
    nt = length(time);
    controls = inputs(:,1:nc);
    wind = controls(:,1);
    wave = controls(:,4);

    states = inputs(:,nc+1:end);

    hf = figure;
    hf.Color = 'w';
    subplot(2,1,1)
    plot(time,wind,'LineWidth',1)
    yline(mean(wind),'--')
    xlabel('Time')
    ylabel('Wind Speed')

   
    subplot(2,1,2)
    plot(time,wave,'LineWidth',1)
    
    xlabel('Time')
    ylabel('Wave Elevation')
    
    x0 = states(1,:);
    
    u_pp = interp1(time,controls,'pchip','pp');
    u_fun = @(t) ppval(u_pp,t);
    
    tic
    [T,X] = ode45(@(t,x) odefun(t,x,A_fun,B_fun,u_fun),time,x0);
    test_time(end+1) = toc;
    
    hf = figure;
    hf.Color = 'w';
    for i = 1:ns/2
    
        subplot(3,1,i)
        hold on;
        plot(time,X(:,i),'Color','#1f77b4','LineWidth',1)
        plot(time,states(:,i),'Color','#ff7f0e','LineWidth',1)
        

        ylabel(reqd_states{i});

    
    
    end
    xlabel('Time')
    legend('DFSM','OpenFAST')
    
    % simulate outputs

    n_out = length(reqd_outputs);
    
    % initialize
    Y = zeros(nt,n_out);

    % loop through and simulate outputs
    for i = 1:nt

        t = time(i);

        u = u_fun(t);
        x = X(i,:);

        w = u(1);

        C_ = squeeze(C_fun(w)); D_ = squeeze(D_fun(w));

        Y(i,:) = C_*x' + D_*u;


    end
    

    hf = figure;
    hf.Color = 'w';
    for i = 1:n_out
    
        subplot(n_out,1,i)
        hold on;
        plot(time,Y(:,i),'Color','#1f77b4','LineWidth',1)
        plot(time,outputs(:,i),'Color','#ff7f0e','LineWidth',1)
        

        ylabel(reqd_outputs{i});

    
    
    end
    xlabel('Time')
    legend('DFSM','OpenFAST')
    
    hf = figure;
    hf.Color = 'w';
    commonFigureProperties

    subplot(2,1,1)
    hold on;
    title('Generator Speed ($\omega_g$) [rad/s]','Interpreter','latex','FontSize',12)
    plot(time,X(:,3),'LineWidth',1.2,'Color','#1f77b4')
    plot(time,states(:,3),'Color','#ff7f0e','LineWidth',1.2)
    taskflag = 'axes'; commonFigureTasks;

    subplot(2,1,2)
    hold on;
    title('Tower-Base Moment ($T_M$) [kNm]','Interpreter','latex','FontSize',12)
    plot(time,Y(:,2),'LineWidth',1.2,'Color','#1f77b4')
    plot(time,outputs(:,2),'Color','#ff7f0e','LineWidth',1.2)
    xlabel('Time [s]','FontSize',12,'Interpreter','latex')
    legend('DFSM','OpenFAST','FontSize',12,'Interpreter','latex','numColumns',2)
    taskflag = 'axes'; commonFigureTasks;

    savefig('DFSM-val.svg')

    %inputs_dfsm = [controls,X];


end

return

function dx = odefun(t,x,A_fun,B_fun,u_fun)

u = u_fun(t);

w = u(1);

A = squeeze(A_fun(w));
B = squeeze(B_fun(w));

dx = A*x + B*u;


end

function sim_details = reorganize_data(sim_details)

    % fid ht eunique load cases
    w_unique = unique(vertcat(sim_details(:).w_mean));
    n_unique = length(w_unique);
    
    % get total load cases
    n_total = length(sim_details);
    
    % rearrange
    if n_unique > 1
        sim_details = reshape(sim_details,[n_total/n_unique,n_unique]);
    end

end
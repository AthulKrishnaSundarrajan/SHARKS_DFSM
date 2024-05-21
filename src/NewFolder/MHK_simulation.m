clc; clear;% close all;

rng(34534)

% file name

root_path = fileparts(which('INSTALL_DFSM'));
data_path = fullfile(root_path,'data');%,'FOWT_rated');
fol_name_cell = {'DFSM_1p6_2'};%{'12','13','14','15','16','17','18','19','20','21','22','23','24'};

nfol = length(fol_name_cell);
AB_cell = cell(nfol,1);
A_cell = cell(nfol,1);

for ifol = 1:nfol

    fol_name = fol_name_cell{ifol};

    sim_path = fullfile(data_path,fol_name);
    
    % file names
    prefix = 'IEA_w_TMD_';
    suffix = '.outb';
    
    % sim_files 
    sim_dir  = dir(fullfile(sim_path,[prefix,'*',suffix]));
    
    nsim = length(sim_dir);%nsim = 10;
    
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
    
    % required states
    reqd_states = {'PtfmPitch','GenSpeed'};
    
    %reqd_states = {'PtfmSurge','PtfmSway','PtfmHeave','PtfmRoll','PtfmPitch','PtfmYaw','GenSpeed'};% ,'GenSpeed','YawBrTAxp'
    
    % filtering arguments for states
    filter_args.filt_states = [~true,~true,false];
    filter_args.filt_states_tf = [1,1,1]; % 2,2 works
    
    % required controls
    reqd_controls = {'RtVAvgxh','GenTq','BldPitch1','Wave1Elev'};
    nc = length(reqd_controls);
    
    % filtering arguments for controls
    filter_args.filt_controls = [~true,false,false];
    filter_args.filt_controls_tf = [1,0,0];
    
    % filter flag
    filter_flag = false;
    
    reqd_outputs = {'TwrBsFxt','TwrBsMyt','NcIMURAys'};
    
    % time
    tmin = 00;
    tmax = [];
    add_dx2 = true;
    
    % extract
    sim_details = simulation_details(sim_files,reqd_states,reqd_controls,reqd_outputs,false,filter_flag,filter_args,add_dx2,tmin,tmax);
    
    split = [0.2,0.8];
    
    % dfsm options
    dfsm_options.ltype = 'LTI';
    dfsm_options.ntype = '';
    dfsm_options.lsamples = nan;
    
    dfsm_options.sampling_type = 'KM';
    dfsm_options.train_test_split = split;
    dfsm_options.scale_flag = 0;
    
    nsamples = 30;
    results_cell = cell(6,1);
    
    
    dfsm_options.nsamples = 10;
    % extract dfsm
    dfsm = DFSM(sim_details,dfsm_options);

    AB = dfsm.deriv.AB;

    for n = 8
        % test index
        ind_test = n;
        plot_flag = true;
        sim_flag = ~false;
        time = sim_details(1).time;
        controls = sim_details(n).controls;
        
        [dfsm_mf,X_cell,dx_cell,Y_cell] = test_dfsm(dfsm,sim_details(ind_test),ind_test,plot_flag,sim_flag);
    
        %results_cell = {time,controls,X_cell,dx_cell,Y_cell};
    end
    
    AB = AB';
    AB_cell{ifol} = AB;
    A = AB(:,nc+1:end);
    A_cell{ifol} = A;


 end

%save('WEIS-P2-Q2.mat','results_cell')

if strcmpi(dfsm_options.ltype,'LPV')

    AB_interp = dfsm.deriv.AB_interp;

    W = dfsm.deriv.W;

    A = AB_interp(:,nc+1:end,:);

    nW = length(W);

    n = 3;

    hf = figure;
    hf.Color = 'w';

    hold on;
    xlabel('Real');ylabel('Imag')
    xline(0)

    ind_stable = false(nW,1);

    for i = 1:nW

        A_ = squeeze(A(i,:,:))';

        eig_A = eig(A_);
        %eig_A = eig_A(n);

        if all(real(eig_A) < 0)
            ind_stable(i) = true;
        end

        plot(real(eig_A),imag(eig_A),'.','MarkerSize',20,'color','r')



    end


end

[nx,nx] = size(A_cell{1});

w = 12:1:24;

train_ind = [1,3,5,7,9,11,13];
test_ind =[2,4,6,8,10,12];



A = zeros(nfol,nx,nx);

for i = 1:nfol
    A(i,:,:) = A_cell{i};

end


A_train = A(train_ind,:,:);
A_test = A(test_ind,:,:);


A_pp = interp1(w(train_ind),A_train,'pchip','pp');
A_op = @(w) ppval(A_pp,w);

n_interp = 500;
w_interp = linspace(12,24,n_interp);

close all;

hf = figure;
hf.Color = 'w';
hold on;
n = 1;

xlabel('Real');ylabel('Imag')

for i = 1:length(train_ind)

    A_ = squeeze(A_train(i,:,:));

    eig_A = eig(A_);
    eig_A = eig_A(n);

    plot(real(eig_A),imag(eig_A),'.','MarkerSize',20,'color','r')


end


for i = 1:length(test_ind)

    A_ = squeeze(A_test(i,:,:));

    eig_A = eig(A_);
    eig_A = eig_A(n);

    plot(real(eig_A),imag(eig_A),'.','MarkerSize',20,'color','k')


end

for i = 1:n_interp 

    w_i = w_interp(i);

    A_i = A_op(w_i);

    eig_A = eig(A_i);
    eig_A = eig_A(n);

    plot(real(eig_A),imag(eig_A),'.','MarkerSize',5,'color','b')



end



AB_r = dfsm.deriv.AB;

filename = 'linear-models.mat';

mat_name = 'DFSM_MHK_R_nowave.mat';

%save(mat_name,'AB_r')

%save(filename,"AB_r",'-append');

return

% extract the linear model
% get the simulation data
[~,~,~,inputs,state_dx,outputs] = sample_data(sim_details(1:4),'KM',5);

AB = linsolve(inputs,state_dx)';
nc = length(reqd_controls);
ns = 2*length(reqd_states);

% initial parameter
par0 = AB(ns/2+1:end,:);

% hard code values that give a stable model
par0_ = [7.41654e-05	-0.00014362	7.6148e-05	-0.00071188	-0.000453151	6.63828e-08	-2.38556	0.00149416	-7.50307e-06	-0.00876847
270.914	-15.2251	-6.20374	15.8834	22.7741	-0.607913	-213.972	-67.391	-0.517846	19.8335
-0.565162	-0.032696	0.0662398	-0.684253	0.043911	0.00231418	-1.2089	1.72458	-0.00496775	-9.64835];


x0 = par0(:);
nx = length(x0);

%%

%
time = sim_details(1).time; dt = time(2) - time(1);
controls = inputs(:,1:nc);
% states = inputs(:,9);


%% testing out GA 
options_hybrid = optimoptions('fmincon','Display','iter','Algorithm','interior-point','MaxIterations',20,'MaxFunctionEvaluations',10000,'UseParallel',true);

% ga options
gaopt = optimoptions('ga','UseParallel',true,'Display','iter','MaxGenerations',nx*2); %,'HybridFcn',{@fmincon,options_hybrid});

% run hybrid optimization scheme
x = ga(@(x) evaluate_loss(x,ns,nc,inputs,state_dx),nx,[],[],[],[],[],[],@(x)test_stability(x,ns,nc),[],gaopt) ; %@(x)findEIG(x,ns,nc)

% get the linear model

%opt_options.Advanced.StepTolerance = 1e-6;
opt_options.TolFun = 1e-6;
opt_options.TolX = 1e-6;
opt_options.MaxIter = 150;
opt_options.Algorithm = 'interior-point';

opt = greyestOptions;
opt.Display ='on';
opt.EnforceStability = true;
opt.SearchMethod ='fmincon';
opt.Focus = 'simulation';
opt.SearchOptions.MaxIterations = 15; 
opt.SearchOptions.FunctionTolerance = 1e-3;
%opt.OutputWeight = diag([1e-0,1e-0,1e-0,1e-0,1e-0,1e-0]);
opt.SearchOptions = opt_options;

% data


% initial system
[A,B,C,D] = LTI_function(x,[],ns,nc);
init_sys = ss(A,B,C,D);

% model
m = idgrey('LTI_function',{'x',x},'c',{ns,nc},0);

data = iddata(state_dx,controls,dt);



options = greyestOptions('Focus','simulation','EnforceStability',true,'Display','on','SearchMethod','fmincon','SearchOptions',opt_options);
% evaluate system
sys = greyest(data,m,opt);

% construct system
dfsm.deriv.AB = [sys.B,sys.A]';



% time and control inputs
time = sim_details(ind_test).time;
controls = sim_details(ind_test).controls;


% test the dfsm
[dfsm_mf,X_cell,dx_cell,Y_cell] = test_dfsm(dfsm,sim_details(ind_test),ind_test,plot_flag,sim_flag);

% extract
dx_m = dx_cell{1}; dx_e = dx_cell{2};
x_m = X_cell{1}; x_e = X_cell{2};

W = opt.OutputWeight;

dx_V = estimate_loss(dx_m,dx_e,W);

dx_training = inputs*dfsm.deriv.AB;

loss_dx = estimate_loss(state_dx,dx_training,eye(ns));

dx_error = state_dx - dx_training;

return






clc; clear; close all;

rng(34534)

% file name

root_path = fileparts(which('INSTALL_DFSM'));
data_path = fullfile(root_path,'data');
fol_name_cell = 'DFSM_1p6_2';



sim_path = fullfile(data_path,fol_name_cell);

% file names
prefix = 'IEA_w_TMD_';
suffix = '.outb';

% sim_files 
sim_dir  = dir(fullfile(sim_path,[prefix,'*',suffix]));

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

% required states
%reqd_states = {'PtfmPitch','GenSpeed','YawBrTAxp'};

reqd_states = {'PtfmPitch','GenSpeed'};% ,'GenSpeed','YawBrTAxp'

% filtering arguments for states
filter_args.filt_states = [~true,~true,false];
filter_args.filt_states_tf = [1,1,1]; % 2,2 works

% required controls
reqd_controls = {'Wind1VelX','GenTq','BldPitch1','Wave1Elev'};
reqd_outputs = {};%'TwrBsMyt','TwrBsFxt','YawBrTAxp','NcIMURAys','GenPwr'};%{'TwrBsFxt','TwrBsMyt','GenPwr','YawBrTAxp'};


% filtering arguments for controls
filter_args.filt_controls = [~true,false,false];
filter_args.filt_controls_tf = [1,0,0];

% filter flag
filter_flag = false;

%reqd_outputs = {};

% time
tmin = 00;
tmax = [];
add_dx2 = true;

% extract
sim_details = simulation_details(sim_files,reqd_states,reqd_controls,reqd_outputs,false,filter_flag,filter_args,add_dx2,tmin,tmax);


split = [0.4,0.6];

% dfsm options
dfsm_options.ltype = 'LTI';
dfsm_options.ntype = '';
dfsm_options.lsamples = nan;

dfsm_options.sampling_type = 'KM';
dfsm_options.train_test_split = split;
dfsm_options.scale_flag = 0;

nsamples = [5];
results_cell = cell(6,1);

[~,~,~,inputs,state_dx,outputs] = sample_data(sim_details(1:8),'KM',5);
outputs = state_dx;

% extract wind speed
wind = inputs(:,1);
wmin = min(wind)+1;
wmax = max(wind)-1;

% get the maximum and minimum values of the wind speed
maxW = max(wind);
minW = min(wind);

% get the size of inputs and outputs
ninputs = size(inputs,2);
noutputs = size(outputs,2);

% create a grid of wind speed values
nW = 3;
W = [11,12,13,14,15,16,17,18,19,20,21,22];nW = length(W);%linspace(wmin,wmax,nW);
data_cell = cell(nW,2);

% initialize offset
offset = 3;

% initialize storage array
IW_ = zeros(length(W),1);

hf = figure;
hf.Color = 'w';
hold on;

hist(wind)
xline(W)

AB_interp = zeros(nW,ninputs,noutputs);

for k = 1:length(W)

    m = W(k);

    if m+offset > maxW
        offset = maxW - m;
    end

    if m-offset < minW
        offset = m-minW;
    end

    % find all values within small range
    Iw = (m+offset >= wind) & (m-offset <= wind);
    IW_(k) = sum(Iw);

    % extract
    input_ = inputs(Iw,:);
    output_ = outputs(Iw,:);

    data_cell{k,1} = input_;
    data_cell{k,2} = output_;

    AB_interp(k,:,:) = linsolve(input_,output_);

end

nc = length(reqd_controls);

A = AB_interp(:,nc+1:end,:);

nW = length(W);

n = 3;

hf = figure;
hf.Color = 'w';

hold on;
xlabel('Real');ylabel('Imag')
xline(0)

ind_stable = false(nW,1);

n = 1;
for i = 1:nW

    A_ = squeeze(A(i,:,:))';

    eig_A = eig(A_);
    

    if all(real(eig_A) < 0)
        ind_stable(i) = true;
    end
    eig_A = eig_A(n);
    plot(real(eig_A),imag(eig_A),'.','MarkerSize',20,'color','r')
    text(real(eig_A)+0.0,imag(eig_A),num2str(W(i)))



end
return
% script to run DFSM NN with PCA

clc; clear; close all

fol_name = 'sim_17ms_600';
t_f = 2;

%% Options for different simulation results
% prefix and suffix of output files
prefix = 'lin_'; % simulations without TTDspFA
suffix = '.outb';


% required states and controls
req_states = {'PtfmPitch','TTDspFA','GenSpeed'}; % 'TTDspFA',
req_controls = {'RtVAvgxh','GenTq','BldPitch1'};

% load
data = loadOFsim(fol_name,prefix,suffix,req_states,req_controls);


% specify which states derivative is included
dindex = 1:length(req_states);

% approximate state derivatives
data = approximateStateDerivatives(data,dindex,1,t_f);

% perform PCA and scale
data = data(1);
data = performPCA(data);

% add state names
state_names = req_states;

for k = 1:length(dindex)
    state_names{end+1} = ['d',req_states{dindex(k)}];
end


X = data.pca_score;
Y = data.state_derivatives;

net = nntrain(X,Y);

U = griddedInterpolant(data.time,data.inputs,'spline');

% simulation options
Time = vertcat(data.time);
t0 = min(Time);
td = max(Time)-min(Time);
TSPAN = [t0 t0+td];
Y0 = data.states(t0 == Time,:);
OPTIONS = odeset('RelTol',1e-5,'AbsTol',1e-5);

% run simulation
[TOUT,YOUT] = ode45(@(t,x) dxfun_NN(t,x,net,U,data.nx,data.nu,data.pca_mu,data.scale_max,data.pca_coeff),TSPAN,Y0,OPTIONS);


 % plot 
hf = figure;
hf.Color = 'w';
hf.Position = [1000 700 560 640];

%sgtitle(title_name)
nl = length(state_names);

for i = 1:nl

subplot(nl,1,i); hold on
plot(TOUT,YOUT(:,i),'r')
plot(data.time,data.states(:,i),'k')
title(state_names{i})
xlim([t0 t0+td])

end
return

function dx = dxfun_NN(t,x,net,U,nx,nu,mu,scale_max,coeff)

% get inputs
u = U(t);

% scale
scale_x = scale_max(1:nx);
scale_u = scale_max(nx+1:nx+nu);

x = x./scale_x';
u = u./scale_u;


X = [x',u];
X = (X-mu)*coeff;

dx = net(X');



end

function data = performPCA(data)

% fitting data
ind = 1;

% extract staes and inputs
states = data(ind).states;
inputs = data(ind).inputs;
nx = size(states,2);data(ind).nx = nx;
nu = size(inputs,2);data(ind).nu = nu;

% combine and scale using maximum value
fitdata = [states,inputs];
scale_max = max(fitdata);
fitdata_scaled = fitdata./scale_max;

% add to data structure
data(ind).fitdata_scaled = fitdata_scaled;
data(ind).scale_max = scale_max;

% perform PCA
[coeff,score,latent,tsquared,explained,mu] = pca(fitdata_scaled);

% store coeff and mu
data(ind).pca_coeff = coeff;
data(ind).pca_mu = mu;
data(ind).pca_score = score;


ind96 = cumsum(explained) < 96;


data(ind).ind96 = ind96;

end

%function net = 
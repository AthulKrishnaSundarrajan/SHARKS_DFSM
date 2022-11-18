clc; clear; close all;


%% folder name
fol_name = 'DFSM_multiwind';
t_f = 2;
reduceflag = 0;
%% Options for different simulation results
% prefix and suffix of output files
prefix = 'lin_'; % simulations without TTDspFA
suffix = '.outb';

%% load OF simulation

% required states and controls
req_states = {'PtfmPitch','TTDspFA','GenSpeed'}; % 'TTDspFA',
req_controls = {'RtVAvgxh','GenTq','BldPitch1'}; % 'GenTq',

% load
data = loadOFsim(fol_name,prefix,suffix,req_states,req_controls);


% specify which states derivative is included
dindex = 1:length(req_states);

% add state names
state_names = req_states;

for k = 1:length(dindex)
    state_names{end+1} = ['d',req_states{dindex(k)}];
end

% approximate state derivatives
data = approximateStateDerivatives(data,dindex,1,t_f);


% perform PCA and scale
data1 = data(1);
data1 = performPCA(data1);

% add state names
state_names = req_states;

for k = 1:length(dindex)
    state_names{end+1} = ['d',req_states{dindex(k)}];
end

% generate model using first simulation
[Ai,Bi] = DFSM_LTI(data1.pca_score,data1.ind96,reduceflag,data1.nx,data1.nu,data1.state_derivatives);

U = griddedInterpolant(data1.time,data1.inputs,'spline');

% simulation options
Time = vertcat(data1.time);
t0 = min(Time);
td = max(Time)-min(Time);
TSPAN = [t0 t0+td];
Y0 = data1.states(t0 == Time,:);
OPTIONS = odeset('RelTol',1e-5,'AbsTol',1e-5);

% run simulation
[TOUT,YOUT] = ode45(@(t,x) dxfun_LTI(t,x,U,data1.nx,data1.nu,data1.pca_mu,data1.scale_max,data1.pca_coeff,Ai,Bi),TSPAN,Y0,OPTIONS);


 % plot 
hf = figure;
hf.Color = 'w';
hf.Position = [1000 700 560 640];

%sgtitle(title_name)
nl = length(state_names);

for i = 1:nl

subplot(nl,1,i); hold on
plot(TOUT,YOUT(:,i),'r')
plot(data1.time,data1.states(:,i),'k')
title(state_names{i})
xlim([t0 t0+td])

end

data2 = data(2); % loadOFsim(fol_name2,prefix,suffix,req_states,req_controls);


% approximate state derivatives
%data2 = approximateStateDerivatives(data2,dindex,1,t_f);


U2 = griddedInterpolant(data2.time,data2.inputs,'spline');

% simulation options
Time = vertcat(data2.time);
t0 = min(Time);
td = max(Time)-min(Time);
TSPAN = [t0 800];
Y0 = data2.states(t0 == Time,:);
OPTIONS = odeset('RelTol',1e-5,'AbsTol',1e-5);

% run simulation
[TOUT2,YOUT2] = ode45(@(t,x) dxfun_LTI(t,x,U,data1.nx,data1.nu,data1.pca_mu,data1.scale_max,data1.pca_coeff,Ai,Bi),TSPAN,Y0,OPTIONS);

% plot 
hf = figure;
hf.Color = 'w';
hf.Position = [1000 700 560 640];

%sgtitle(title_name)
nl = length(state_names);

for i = 1:nl

subplot(nl,1,i); hold on
plot(TOUT2,YOUT2(:,i),'r')
plot(data2.time,data2.states(:,i),'k')
title(state_names{i})
xlim([t0 t0+td])

end



return

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

function [Ai,Bi] = DFSM_LTI(score,ind96,reduceflag,nx,nu,state_derivatives)


% transform data

if reduceflag
    fitdata_ = score(:,ind96);
else
    fitdata_ = score;

end

% find the best linear model
ABD = linsolve(fitdata_,state_derivatives);

% transpose
ABD = ABD';

% extract
Ai = ABD(:,1:nx);
Bi = ABD(:,nx+1:nx+nu);


end

function dx = dxfun_LTI(t,x,U,nx,nu,mu,scale_max,coeff,Ai,Bi)

% get inputs
u = U(t);

% scale
scale_x = scale_max(1:nx);
scale_u = scale_max(nx+1:nx+nu);

x = x./scale_x';
u = u./scale_u;


X = [x',u];
X = (X-mu)*coeff;

dx = Ai*(X(1:nx))' + Bi*(X(nx+1:nx+nu))';

end
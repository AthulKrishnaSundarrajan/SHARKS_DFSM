% script to run DFSM LTI with PCA

clc; clear; close all

fol_name = 'sim_17ms_600';
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
dx = data.state_derivatives;
Time = data.time;

dx_pp = pchip(Time,dx');
dx_fun = @(t)ppval(dx_pp,t);


% simulation options

t0 = min(Time);
td = max(Time)-min(Time);
TSPAN = [t0 t0+td];
Y0 = data.states(t0 == Time,:);
OPTIONS = odeset('RelTol',1e-5,'AbsTol',1e-5);

[T,Y] = ode45(@(t,x) odefun(t,x,dx_fun),TSPAN,Y0,OPTIONS);



 % plot 
hf = figure;
hf.Color = 'w';
hf.Position = [1000 700 560 640];

%sgtitle(title_name)
nl = length(state_names);

for i = 1:nl

subplot(nl,1,i); hold on
plot(T,Y(:,i),'r--')
plot(data.time,data.states(:,i),'k')
title(state_names{i})
xlim([t0 t0+td])

end


return

function dx = odefun(t,x,dx_fun)

    dx = dx_fun(t);


end
%
function y = smoothData2(t,x,varargin)

% time step (assumed constant in t)
dt = t(2) - t(1);

% calculate sample rate
h = 1/dt;

% create lowpass IIR response
D = designfilt('lowpassiir', 'FilterOrder', 7, ...
     'PassbandFrequency', 5, 'PassbandRipple', 0.001,...
     'SampleRate', h);

% filter the data
y = filtfilt(D,x);

end

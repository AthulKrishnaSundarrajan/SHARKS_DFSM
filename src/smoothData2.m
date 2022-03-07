function y = smoothData2(t,x,varargin)

dt      = t(2) - t(1);
h = 1/dt;

D = designfilt('lowpassiir', 'FilterOrder', 7, ...
     'PassbandFrequency', 5, 'PassbandRipple', 0.001,...
     'SampleRate', h);

y = filtfilt(D,x);
 

end

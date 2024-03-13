function[fs,FFT_sig] = perform_FFT(time,signal)


nt = length(time);

dtDelta0 = time(2)-time(1);

dt = (max(time) - min(time))/(nt-1);

if dt ~= dtDelta0
    disp('[WARN] dt from tmax-tmin different from dt from t2-t1')
end

Fs = 1/dt;

nFFTAll = 2^ceil(log(nt)*0.99999999999/log(2));

nExp = log(nFFTAll)/log(2) - 1;

nPerSeg = 2^nExp;

window = hamming(nPerSeg,'symmetric');

[FFT_sig,fs] = pwelch(signal,window,[],[],Fs);


end
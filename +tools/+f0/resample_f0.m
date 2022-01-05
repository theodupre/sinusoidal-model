function [tn, fn] = resample_f0(regime, T0, Tend, sampleStart, sampleEnd, fs)

    f = regime(sampleStart:sampleEnd,3)*10000/120;
    t = linspace(T0,Tend,sampleEnd-sampleStart+1)';
    tn = linspace(T0,Tend,fs/1000*(Tend-T0))';
    fn = interp1(t,f,tn);

end
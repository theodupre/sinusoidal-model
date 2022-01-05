clear; close all;

addpath(genpath('./'));

%files = dir('D:\*\processed\2008\**\[Audio]*.wav');
files = dir('./audio/[Audio]test*.wav');

%saveFolder = 'D:\ManipMesures3D\ValidationSon3D\stimuli\';
saveFolder = '.\audio\';

for k = 1:1%length(files)
    
    fprintf('Processing file %s (%i/%i) ... \n', files(k).name, k, length(files))
%     tic
    filek = [files(k).folder '\' files(k).name];
    
    info = audioinfo(filek);
    
%     Tend = floor(info.Duration);
%     if Tend > 240
%         Tend = 240;
%     end
    
    [regime,fs] = audioread(strrep(filek, 'Audio', 'CAN'));

    T0 = 0;
    Tend = 3;
    sampleStart = T0*fs+1;
    sampleEnd = Tend*fs;

    [tn, fn] = tools.f0.resample_f0(regime, T0, Tend, sampleStart, sampleEnd, fs);

    [s, fs] = audioread(filek, [sampleStart sampleEnd]);

    mcSinusoidal = zeros(sampleEnd - sampleStart+1, info.NumChannels);
    mcResidual = mcSinusoidal;

    rangeThreshold = 5;
    winflag = 6; % Hann:3  Hamming:7  Blackmann:5  Blackmann-Harris:6 

    processFrame = 2;
    partialDisp = false;
    for i = 1:1 %info.NumChannels % Process only one channel for an example
        fprintf('Channel %i/%i \n',i,info.NumChannels)
        for j = 1:(Tend-T0)/processFrame
         
            samples = (1:processFrame*fs)+(j-1)*processFrame*fs;
            [mcSinusoidal(samples,i), mcResidual(samples,i)] = analysisSynthesis(s(samples,i), fn, fs, winflag, rangeThreshold, partialDisp);
            partialDisp = false;
        end
        if samples(end)+ sampleStart < sampleEnd
            samples = samples(end):(sampleEnd-sampleStart);
            [mcSinusoidal(samples,i), mcResidual(samples,i)] = analysisSynthesis(s(samples,i), fn, fs, winflag, rangeThreshold, partialDisp);
        end
    end
    %fprintf('Last Sample Processed %i/%i \n', samples(end), sampleEnd-sampleStart);
    [~,name,~] = fileparts(filek);

    audiowrite([saveFolder, name, '_engine.wav'], mcSinusoidal, fs)
    audiowrite([saveFolder, name, '_noise.wav'], mcResidual, fs)

    %fprintf('Time to process file %i : %f seconds. \n', k, toc);
    
end
function [sinusoidal, residual] = analysisSynthesis(wav, f0, fs, winflag, rangethres, partialDisp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SINUSOIDAL ANALYSIS PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Peak selection
peakselflag = true;

% Track duration
trackdurflag = false;

% Display name of analysis window in the terminal
%fprintf(1,'%s analysis window\n',tools.dsp.infowin(winflag,'name'));

% Flag for causality of first window
causalflag = {'causal','non','anti'};
cf = 3;

% Flag for log magnitude spectrum
logflag = {'dbr','dbp','nep','oct','bel'};
lmsf = 2;

% Flag for parameter estimation
paramestflag = {'nne','lin','log','pow'};
pef = 4;

% Partial tracking flag
ptrackflag = {'','p2p'};
ptf = 2;

% Number of fundamental periods
if peakselflag
    nT0 = 6;
    % Oversampling factor
    osfac = 4;
else
    nT0 = 4;
    osfac = 2;
end

% Normalize analysis window
normflag = true;

% Use zero phase window
zphflag = true;

% Maximum number of peaks to retrieve from analysis
maxnpeak = 200;

% Return MAXNPEAK frequency bins
npeakflag = true;

% Replace -Inf in spectrogram
nanflag = false;

% Peak shape threshold (normalized)
shapethres = 0.8;

% Peak range threshold (dB power)
%rangethres = 15;

% Relative threshold (dB power)
relthres = -90;

% Absolute threshold (dB power)
absthres = -100;

% Minimum partial track segment duration (ms)
durthres = 10.5;

% Connect over (ms)
gapthres = 22.5;

% Resynthesis flag
synthflag = {'OLA','PI','PRFI'};
sf = 2;

% Display resynthesis info
dispflag = false;

% Reference f0
ref0 = tools.f0.reference_f0(f0);

% Frame size = n*T0
framelen = tools.dsp.framesize(f0,fs,nT0);

% 75% overlap
hop = tools.dsp.hopsize(framelen,0.75);

% FFT size
nfft = tools.dsp.fftsize(framelen,osfac);

% Frequency difference for peak matching (Hz)
freqdiff = tools.dsp.freq_diff4peak_matching(ref0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SINUSOIDAL ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[amplitude,frequency,phase,center_frame,npartial,nframe,nchannel,nsample,dc] = sinusoidal_analysis(wav,framelen,hop,nfft,fs,...
    winflag,maxnpeak,shapethres,rangethres,relthres,absthres,durthres,gapthres,freqdiff,...
    causalflag{cf},paramestflag{pef},ptrackflag{ptf},normflag,zphflag,npeakflag,peakselflag,trackdurflag);

if partialDisp
    
    % Figure layout
    fig_layout.font = 'Times New Roman';
    fig_layout.axesfs = 14;
    fig_layout.titlefs = 22;
    fig_layout.bckgdc = [1 1 1];
    fig_layout.cmap = 'jet';
    fig_layout.figsize = [15 10];
    fig_layout.figpos = [0.5 0.5 fig_layout.figsize-0.5];
    fig_layout.figunit = 'centimeters';
    fig_layout.msize = 7;
    fig_layout.mtype = '.';
    fig_layout.linsty = '-'; %'none'
    fig_layout.meshsty = 'row'; %'both'
    fig_layout.disp = 'on';
    fig_layout.print = 'opengl';
    
    % Axes label
    axes_lbl.tlbl = 'Time (s)';
    axes_lbl.flbl = 'Frequency (kHz)';
    axes_lbl.dblbl = 'Spectral Energy (dB)';
    
    % Axes limits
    axes_lim.flim = [0 8];
    axes_lim.dblim = [-120 0];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT PARTIAL TRACKING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Log Mag Amp peaks
    plot_part.specpeak = tools.math.lin2log(amplitude,logflag{lmsf},nanflag);
    
    % Time peaks
    plot_part.time = repmat(center_frame/fs,[1,npartial])';
    
    % Frequency peaks
    plot_part.frequency = frequency/1000;
    
    % Time limits
    axes_lim.tlim = [plot_part.time(1,1) plot_part.time(1,end)];
    
    % Title
    axes_lbl.ttl = 'Partial Tracking';
    
    % Make figure
    tools.plot.mkfigpeakgram(plot_part,axes_lim,axes_lbl,fig_layout);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SINUSOIDAL RESYNTHESIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sinusoidal,partial,amp_partial,freq_partial,phase_partial] = sinusoidal_resynthesis(amplitude,frequency,phase,...
    framelen,hop,fs,nsample,center_frame,npartial,nframe,nchannel,durthres,gapthres,freqdiff,...
    winflag,causalflag{cf},synthflag{sf},ptrackflag{ptf},~trackdurflag,dispflag);

% Make residual
residual = wav - sinusoidal;

% Normalize residual to -16dB
%normres = tools.wav.normdb(residual,-16,'rms');

% Calculate signal-to-resynthesis energy ratio (SRER)
% srer_db = tools.wav.srer(wav,residual);
% 
% % Time in signal reference (s)
% time = tools.plot.mktime(nsample,fs);


end
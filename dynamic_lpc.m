function [A,e,center_frame,short_fft,lpc_env] = dynamic_lpc(wav,hop,framelen,winflag,causalflag,lpc_order)
%DYNAMIC_LPC Frame-by-frame linear prediction coefficients.
%
%   [A,E,CFRAME,FFT_FR,LPC_ENV] = DYNAMIC_LPC(WAV,H,M,WINFLAG,CENTERWIN,ORDER)

% 2016 MCaetano (Revised)
% 2019 MCaetano SMT 0.1.0
% 2020 MCaetano SMT 0.1.1 (Revised)
% 2020 MCaetano SMT 0.2.0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WORKING PROPERLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LPC frame by frame

nfft = 2^nextpow2(2*framelen);

normflag = false;

[time_frame,center_frame,nsample,nframe,nchannel,dc] = sof(wav,hop,framelen,winflag,causalflag,normflag);

short_fft = fft(time_frame,nfft);

short_fft(short_fft==0) = realmin;

% Raw autocorrelation coefficients: deterministic signal
R = ifft(abs(short_fft).^2);

R(R==0) = realmin;

% Biased autocorrelation estimate: Random (stochastic) process
% R = R./framelen;

[A,e,k] = levinson(R,lpc_order);

A(isnan(A)) = realmin;

e(isnan(e)) = realmin;

% Output LPC in columns
A = A';

% calculate spectral envelope curve from linear prediction coefficients
aux_lpc = abs(1./fft(A,nfft));

lpc_env = sqrt(e').*aux_lpc(1:nfft/2+1,:);

end

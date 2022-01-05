function [filtered,filt_sig] = dynamic_filter(exct_sig,lpc_coeffs,lpc_err,frame_pos,flag)
%
% LPC estimation frame by frame
%
%   ECT_SIG: excitation (signal to be filtered/estimated)
%   LPC_COEFFS: Filter coefficients estimated by linear prediction
%   (levinson)
%   LPC_ERR: Linear prediction estimation error (squared) or signal energy
%   inside a frame
%   FRAME_POS: Center of analysis/sintesis window
%   FLAG: 1 (one) when EXCT_SIG is white noise and 0 (zero) when EXCT_SIG is
%   original signal
%
%   When EXCT_SIG is original signal (where LPC coefficients come from),
%   FILTERED is whitened EXCT_SIG (inverse
%   filtered by A(z))
%
%   When EXCT_SIG is white noise, FILTERED is white noise
%   colored by filter A(z) (LPC)

% 2016 MCaetano (Revised)
% 2019 MCaetano SMT 0.1.0
% 2020 MCaetano SMT 0.1.1 (Revised)
% 2020 MCaetano SMT 0.2.0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WORKING PROPERLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Filtering:
%     y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                  - a(2)*y(n-1) - ... - a(na+1)*y(n-na)

dur = length(exct_sig);

[lpc_order] = size(lpc_coeffs,1);

lpc_order = lpc_order - 1;

filt_sig = zeros(dur,1); % The filtered waveform

for count = 1:dur
    
    [v,p] = min(abs(count - frame_pos));
    
    
    if count == 1
        
        filt_sig(count) = exct_sig(count);
        
    else
        
        order = min(lpc_order, count-1);
        
        if flag
            
            %   Notice that FILT_SIG is the prediction error because it
            %   contains 1-sum(a_k*s(n-k)) = A(z)
            %
            %   Filtering:
            %     y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
            %                  - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
            %
            %   Below:
            %
            %     y(n) = b(1)*x(n) - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
            
            filt_sig(count) = exct_sig(count) - sum(lpc_coeffs(2:order+1,p).*filt_sig(count-1:-1:count-order));
            
        else
            
            %   The signal estimated by linear prediction is simply -sum(a_k*s(n-k))
            %   Filtering:
            %     y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
            %                  - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
            %
            %
            %   Below:
            %     y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
            
            filt_sig(count) = - sum(lpc_coeffs(2:order+1,p).*exct_sig(count-1:-1:count-order));
            
        end
        
    end
    
end

if flag == 1
    
    % Apply energy for Random (stochastic) process
    ener = exp(interp1(frame_pos, log(sqrt(lpc_err)),[1:dur], 'linear','extrap'))';
    %ener = resample(sqrt(lpc_err),dur,num_frames);
    filtered = ener.*filt_sig;
    
elseif flag == 2
    
    filtered = filt_sig;
    
else
    
    %   Output prediction error (whitened input signal)
    
    %     filtered = filt_sig; %  only outputs linear prediction
    filtered = exct_sig - filt_sig;
    
end

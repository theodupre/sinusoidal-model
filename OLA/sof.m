function [time_frame,nsample,dc,centerframe] = sof(wav,framelen,hop,winflag,cfwflag,normflag)
%SOF Split into overlapping frames.
%   FR = SOF(S,M,H,WINFLAG,CFWFLAG,NORMFLAG) splits the input S into
%   overlapping frames FR of length M with a hop size H.
%
%   WINFLAG is a numeric flag that specifies the following windows
%
%   1 - Rectangular
%   2 - Bartlett
%   3 - Hann
%   4 - Hanning
%   5 - Blackman
%   6 - Blackman-Harris
%   7 - Hamming
%
%   CFWFLAG is a string that determines the center of the first analysis
%   window. CFWFLAG can be 'ONE', 'HALF', or 'NHALF'. The sample CFW
%   corresponding to the center of the first window is obtained as
%   CFW = cfw(M,CFWFLAG).
%
%   NORMFLAG is a boolean that specifies if the window is normalized as
%   normw(n)=w(n)/sum(w(n)). NORMFLAG=TRUE normalizes and FALSE does not.
%
%   [FR,NSAMPLE,DC,CFR] = SOF(...) also returns the original length NSAMPLE
%   of S (in samples), the dc value DC of S, and the vector CFR with the
%   samples corresponding to the center of the time frames FR.
%
%   See also OLA

% 2016 M Caetano
% 2019 MCaetano (Revised)
% 2020 MCaetano SMT 0.1.1 (Revised)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK FUNCTION CALL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check number of input arguments
narginchk(6,6);

% Check number of output arguments
nargoutchk(0,4);

%TODO: CHECK INPUTS (CLASS/VALUE/NAN)
%TODO: HANDLE STEREO SOUNDS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION BODY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make SIG column vector
wav = wav(:);% WARNING! THIS LINE CONCATENATES THE CHANNELS OF STEREO SOUNDS

% Get NSAMPLE
nsample = size(wav,1);

% Prepare SOF
[analysis_window,dc,nframe,centerframe] = sofprep(nsample,framelen,hop,cfwflag,winflag,normflag);

% Execute SOF
[time_frame] = sofexe(wav,nsample,analysis_window,framelen,centerframe,nframe);

end

% FUNCTION THAT PREPARES SOF
function [analysis_window,dc,nframe,centerframe] = sofprep(nsample,framelen,hop,cfwflag,winflag,normflag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK WINDOW OVERLAP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of overlapping samples
noverlap = framelen - hop;

%hop == framelen
if noverlap == 0
    
    warning('SMT:AdjacentFrames',['Frames are adjancent.\n'...
        'Frames do not overlap.\n Typically, consecutive time_frame overlap by 50%.']);
    
    %hop > framelen
elseif noverlap < 0
    
    warning('SMT:NonOverlappingFrames',['Frames do not overlap.\n'...
        'There is a gap of %d samples between consecutive time_frame.'],...
        abs(noverlap));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE REMAINING VARIABLES TO GENERATE FRAMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make analysis window
analysis_window = mkcolawin(framelen,winflag);

% Normalize analysis window energy
if normflag
    % sum(ANALYSIS_WINDOW) = 1
    dc = sum(analysis_window);
    analysis_window = analysis_window/dc;
else
    % No normalization
    dc = 1;
end

% Number of time_frame
nframe = numframe(nsample,framelen,hop,cfwflag);

% Center of each frame in signal reference
centerframe = f2s(1:nframe,cfw(framelen,cfwflag),hop);

end

% FUNCTION THAT EXECUTES SOF
function [windowed_frames,nsample,centerframe] = sofexe(wav,nsample,analysis_window,framelen,centerframe,nframe)

% Initialize the matrix FRAMES that will hold the signal time_frame
time_frame = zeros(framelen,nframe);

% Initialize WINDOWED_FRAMES
windowed_frames = zeros(framelen,nframe);

% For each frame
for iframe = 1:nframe
    
    % Make each frame
    time_frame(:,iframe) = mkframe(wav,centerframe(iframe),framelen,nsample);
    
    % Perform windowing
    windowed_frames(:,iframe) = analysis_window.*time_frame(:,iframe);
    
end

end

% FUNCTION THAT MAKES EACH FRAME
function frame = mkframe(wav,cf,framelen,nsample)

% Initialize frame
frame = zeros(framelen,1);

% First half window (right of CW)
fhw = rhw(framelen);

% Second half window (left of CW)
shw = lhw(framelen);

if cf < 1
    
    up_bound = cf + fhw;
    
    if up_bound >= 1
        
        frame(framelen-(up_bound-1):framelen) = wav(1:up_bound);
        
    end
    
elseif cf > nsample
    
    low_bound = cf - shw;
    
    if low_bound <= nsample
        
        frame(1:nsample-(low_bound-1)) = wav(low_bound:nsample);
        
    end
    
else
    
    % Sample position of lower window bound in signal reference
    low_bound = max(1,cf-shw);
    
    % Sample position of upper window bound in signal reference
    up_bound = min(nsample,cf+fhw);
    
    if cf - shw < 1
        
        frame(framelen-(up_bound-low_bound):framelen) = wav(low_bound:up_bound);
        
    elseif cf + fhw > nsample
        
        frame(1:(up_bound-low_bound)+1) = wav(low_bound:up_bound);
        
    else
        
        frame(1:framelen) = wav(low_bound:up_bound);
        
    end
    
end

end

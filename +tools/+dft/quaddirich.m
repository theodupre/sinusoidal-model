function qdirich = quaddirich(bin,framelen,nfft,coeff0,coeff1,coeff2,coeff3,normflag,zphflag)
%QUADDIRICH Triple dirichlet kernel.
%   DD = QUADDIRICH(BIN,FRAMELEN,NFFT,C1,C2,C3,NORMFLAG,ZPHFLAG) returns the
%   quadruple dirichlet kernel DD = C1*D(BIN) + C2*[D(BIN-S/2) + D(BIN+S/2)],
%   where D(BIN) is the Dirichlet kernel of degree FRAMELEN evaluated at
%   BIN and sampled at NFFT generated by TOOLS.DFT.DIRICHLET. C1 and C2 are
%   the window-dependent coefficients applied to the center Dirichlet
%   lernel and shifted Dirichlet kernels respectively. S is the shift given
%   by S = 2*NFFT/FRAMELEN. NORMFLAG is a ligocal flag that determines if D
%   is normalized by FRAMELEN. NORMFLAG = TRUE sets normalization and
%   NORMFLAG = FALSE does not. ZPHFLAG is a logical flag that determines if
%   D is zero phase or linear phase. ZPHFLAG = TRUE sets zero-phase and
%   ZPHFLAG = FALSE sets linear-phase.
%
%   See also DIRICHLET

% 2021 M Caetano SMT
% $Id 2021 M Caetano SM 0.9.0-alpha.1 $Id


% TODO: CHECK ATTRIBUTES OF INPUT ARGUMENTS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check the number of input arguments
narginchk(9,9);

% Check the number of output arguments
nargoutchk(0,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Multiplicative signal for causal Hann window
if zphflag
    
    sig = 1;
    
else
    
    sig = -1;
    
end

% SHIFT==Position of first zero-crossing of rectangular window
if tools.misc.isodd(framelen)
    
    % Frequency bin shift
    shift = tools.spec.interp_bin(nfft,framelen-1);
    
    % Normalization factor
    normfac = coeff0*framelen - coeff1 + coeff2 - coeff3;
    
else
    
    % Frequency bin shift
    shift = tools.spec.interp_bin(nfft,framelen);
    
    % Normalization factor
    normfac = coeff0*framelen;
    
end

% Non-shifted Dirichlet kernel
dirich_center = tools.dft.dirichlet(bin,framelen,nfft,false,zphflag);

% Dirichlet kernel shifted left by one frequency bin
dirich_left = tools.dft.dirichlet(bin + shift,framelen,nfft,false,zphflag);

% Dirichlet kernel shifted left by two frequency bins
dirich_left2 = tools.dft.dirichlet(bin + 2*shift,framelen,nfft,false,zphflag);

% Dirichlet kernel shifted left by three frequency bins
dirich_left3 = tools.dft.dirichlet(bin + 3*shift,framelen,nfft,false,zphflag);

% Dirichlet kernel shifted left by one frequency bin
dirich_right = tools.dft.dirichlet(bin - shift,framelen,nfft,false,zphflag);

% Dirichlet kernel shifted left by two frequency bins
dirich_right2 = tools.dft.dirichlet(bin - 2*shift,framelen,nfft,false,zphflag);

% Dirichlet kernel shifted left by three frequency bins
dirich_right3 = tools.dft.dirichlet(bin - 3*shift,framelen,nfft,false,zphflag);

% Final quadruple Dirichlet kernel
qdirich = coeff0*dirich_center + sig*coeff1*0.5*(dirich_right + dirich_left) + ...
    coeff2*0.5*(dirich_right2 + dirich_left2) + sig*coeff3*0.5*(dirich_right3 + dirich_left3);

if normflag
    
    qdirich = qdirich/normfac;
    
end

end

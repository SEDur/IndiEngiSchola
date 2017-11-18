function [alpha,one_minus_alpha]=tau2alpha(tau_ms,fs,bitdepth)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% translate time constant to alpha...
% INPUTS:   tau (in mSec) 
%           fs  sample rate (Hz)
%           bitdepth [optional] in bits. for a 1.23bit number type 23 ...
% OUTPUTS:  alpha  (decayfactor 0:1)
%           one_minus_alpha (amplitudefactor 1:0)
%
% Sean Thomson  B&W Group Ltd, Jan 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fc=1/(2*pi*(tau_ms/1000));
alpha=exp(-2*pi*fc/fs); %decayfactor
one_minus_alpha=1-alpha; %amplitudefactor

if nargin == 3
        alpha=          round(          alpha*2^(bitdepth)) * 2^-(bitdepth);
        one_minus_alpha=round(one_minus_alpha*2^(bitdepth)) * 2^-(bitdepth);
end

if nargin == 3
if alpha==1 || alpha==0 || one_minus_alpha==1 || one_minus_alpha==0
    disp(['!!! rounding limit reached:  alpha= ',num2str(alpha),' , one_minus_alpha= ',num2str(one_minus_alpha)]);
    error('alpha precision limits!!')
end
end

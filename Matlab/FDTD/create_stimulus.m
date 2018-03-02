function [stimulus, infoString] = create_stimulus(strType, A, fc, T, numZeros, LPFfreq, fs)
% function [stimulus, infoString] = create_stimulus(strType, A, fc, T, numZeros, LPFfreq, fs)
%
% DESCRIPTION:
% Creates a time domain signal normally used as a stimulus input to a simulated 
% dynamic eq system. All args required but numZeros, LPFfreq can be empty arrays
%
% INPUTS: 
%   strType:  one of these strings 
%        'tone'  : constant tone of freq "fc" and amplitude "A" 
%        'burst' : tone of multi freq "fc" and amplitue "A" spread among "T"
%        'ramp'  : tone freq "fc", ramp up to A down to 0,   
%        'sweep' : sweep from fc(1) to fc(2) fc=[f1 f2]
%        'noise' : randn of gain A
%        'tonedc': tone Fc gain A(1) plus DC of amplitude A(2)
%   A:  amplitude
%   fc: target frequency in Hz (of the tone, ramp etc.) 
%   T:  signal duration in seconds 
%   numZeros: [optional pass empty] number ofzeros to prepend
%   LPFfreq: [optional pass empty] filter the signal with a 2nd order butter of LPFfreq in Hz
%   fs: sample rate in Hz
% 
% OUTPUTS:
%   stimulus: to signal
%   infoString: some information describing the signal
%
% TODO: this could live in matlabfiles, and be re-writen
% Sean Thomson, Aug 2017    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

switch strType
    case 'tone'  %  constant tone of fc and A 
        %signal generation
        signalLength = round(T * fs);
        t = [0:signalLength-1]' / fs;      % array of time
        stimulus = A * sin(2 * pi * fc .* t);
        infoString = ['StimuTone_', num2str(fc), 'Hz_A', num2str(A,2), '_T', num2str(T), 's'];
    case 'burst' %  tone of multi fc and A spread among T
        signalLength = round(T * fs);
        t = [0:signalLength-1]' / fs;      % array of time
        cols = length(A);
        assert(0 == rem(signalLength, cols),'we need 0==rem(N,cols) ')
        tm = reshape(t, signalLength / cols, cols);
        stimulus = [];
        for n=1 : cols
            if isscalar(fc)                
                x = A(n) * sin(2 * pi * fc .* tm(:,n)); %constant fc
            else
                x = A(n) * sin(2 * pi * fc(n) .* tm(:,n));%changing fc
            end                
            stimulus = [stimulus; x];
        end
        infoString = ['StimuBurst_', num2str(fc), 'Hz_A', num2str(A,2), '_T' ,num2str(T),'s'];
    case 'ramp'  %  ramp up down          
        signalLength = round(T * fs);
        t = [0 : signalLength-1]' / fs;      % array of time
        As = [linspace(0, A, floor(signalLength / 2))';
            linspace(A, 0, ceil(signalLength / 2))'];
        stimulus = As .* sin(2 * pi * fc .* t);
        infoString = ['StimuRamp_', num2str(fc), 'Hz_A', num2str(A, 2), '_T', num2str(T), 's'];
    case 'sweep' %  sweep    fc=[f1 f2]
        stimulus = sweep_create(fc(1), fc(2), fs, T, 2*numZeros / fs, A); 
        infoString = ['StimuSweep_', num2str(fc), 'Hz_A', num2str(A, 2), '_T', num2str(T), 's'];
    case 'noise' % noise (LP filtered to Fc Hz)
        signalLength = round(T * fs);
        t = [0 : signalLength-1]' / fs;      % array of time
        %rng(0);
        stimulus = A * randn(size(t));
        infoString = ['StimuNoise_A', num2str(A), '_T', num2str(T), 's'];
    case 'tonedc'  %  constant tone of fc and A 
        %signal generation
        signalLength = round(T * fs);
        t = [0 : signalLength-1]' / fs;      % array of time
        stimulus = A(1) * sin(2 * pi * fc .* t) + A(2);
        infoString = ['StimuToneDC_', num2str(fc), 'Hz_A', num2str(A(1), 2),...
            'plus', num2str(A(2), 2), 'DC_T', num2str(T), 's'];        
    otherwise 
        error('unknown stimulus type')
        % TODO: Add signal with unitary amplitude and random phase
end

% band limit because it's a woofer
if isfinite(LPFfreq)
    [b,a] = butter(2, LPFfreq / (fs / 2));
    stimulus = filter(b, a, stimulus);
end

%pre pend zeros
stimulus = [zeros(numZeros, 1); stimulus(1 : end - numZeros)]; 



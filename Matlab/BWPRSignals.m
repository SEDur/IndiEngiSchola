% BW Noise PSD Prototype Script
%
% SD, ST, JV
%
% Last Access
% 5/2017

classdef BWPRSignals
    properties
        fftLength
        trials
        samplesPerFrame
        randSeed
        
    end
    properties (SetAccess = private)
        NoisePowerSpectralDensity
        fs
        arithmetic
        iwl
        ifl
        currentFrame
        totalFrames
    end
    properties (Dependent)
        signal
    end
    methods
        function obj = BWPRSignals (varargin)
            % BWNoisePSD Object constructor
            if nargin == 1
                % Instantiate object properties from struct
                % Setup Noise PSD Params
                obj.fftLength = 2048;
                obj.trials = 10;
                obj.randSeed = 1;
                obj.fs = 96000;
                obj.frequencyScale = 'fs';
                obj.arithmetic = 'double';
                obj.iwl = 8;
                obj.ifl = 25;
                
            else
                % Declare Parser
                p = inputParser;
                % Def default property values
                % Setup Noise PSD Params
                defaultfftLength = 2048;
                defaulttrials = 10;
                defaultrandSeed = 1;
                defaultfs = 96000;
                defaultfrequencyScale = 'fs';
                defaultarithmetic = 'double';
                defaultiwl = 8;
                defaultifl = 24;
                
                % Def Func Inputs
                % Setup Noise PSD Params
                addOptional(p, 'fftlength', defaultfftLength, @isnumeric);
                addOptional(p, 'trials', defaulttrials, @isnumeric);
                addOptional(p, 'randSeed', defaultrandSeed, @isnumeric);
                addOptional(p, 'fs', defaultfs, @isnumeric);
                addOptional(p, 'frequencyscale', defaultfrequencyScale);
                addOptional(p, 'arithmetic', defaultarithmetic);
                addOptional(p, 'iwl', defaultiwl);
                addOptional(p, 'ifl', defaultifl);
                
                % Parse Inputs
                parse(p,varargin{:});
                % Apply parsed inputs to instantiate object properties
                obj.fftLength = p.Results.fftlength;
                obj.fs = p.Results.fs;
                obj.randSeed = p.Results.randSeed;
                obj.trials = p.Results.trials;
                obj.frequencyScale = p.Results.frequencyscale;
                obj.arithmetic = p.Results.arithmetic;
                obj.iwl = p.Results.iwl;
                obj.ifl = p.Results.ifl;
            end
            
            [stimTimeDomain, stimFreqDomain] = obj.makeStimulus();
            filterOutFreqDomain = obj.getFilterResponse(stimTimeDomain, myFilter);
            obj.NoisePowerSpectralDensity = obj.theMagic(filterOutFreqDomain, stimFreqDomain);
            release(myFilter);
            
            % Scale Pnn to be a true PSD
            switch obj.frequencyScale
                case 'norm'
                    obj.NoisePowerSpectralDensity = obj.NoisePowerSpectralDensity/(2*pi);
                case 'fs'
                    obj.NoisePowerSpectralDensity = obj.NoisePowerSpectralDensity/obj.fs;
                otherwise
                    obj.NoisePowerSpectralDensity = obj.NoisePowerSpectralDensity/obj.fs;
            end
            
            % Convert a 'twosided' spectrum (and frequencies) to a 'onesided' spectrum.
            if rem(ceil(obj.fftLength/2),2)
                select = 1:(obj.fftLength+1)/2;                    % ODD;  take only [0,pi)
                obj.NoisePowerSpectralDensity = [obj.NoisePowerSpectralDensity(1,:); 2*obj.NoisePowerSpectralDensity(select(2:end),:)]; % Don't multiply DC term by 2.
            else
                select = 1:obj.fftLength/2+1;                      % EVEN; take only [0,pi]
                obj.NoisePowerSpectralDensity = [obj.NoisePowerSpectralDensity(1,:); 2*obj.NoisePowerSpectralDensity(select(2:end-1),:); obj.NoisePowerSpectralDensity(select(end),:)]; % Don't multiple DC & Nyquist by 2.
            end
            obj.plot(obj);
        end
        function value = get.totalRoundOffNoisePower(obj)
            range = [0 obj.fs/2];
            value = bandpower(obj.NoisePowerSpectralDensity,linspace(0,48000,ceil(obj.fftLength/2)+1),range,'psd');
        end
    end
    methods (Static)
        function plot(obj)
            semilogx(linspace(0,obj.fs,ceil(obj.fftLength/2)+1),db(obj.NoisePowerSpectralDensity,'p'));
            grid;
            title('round off noise psd')
            ylabel('[dB/Hz]')
            display(strcat('The total RoundOff Noise Power is_',num2str(obj.totalRoundOffNoisePower),' dB'))
        end
    end
    methods (Access = protected)
        function [stimTimeDomain, stimFreqDomain] = makeStimulus(obj)
            % Computer Noise PSD Spectrum
            
            rng(obj.randSeed);
            % Generate Stimulus
            if rem(obj.fftLength,2) == 0 % If even length
                periodicSignal = 2*pi*rand(obj.fftLength/2-1,obj.trials);
                paddedPeriodicSignal = [zeros(1,obj.trials); periodicSignal;...
                    zeros(1,obj.trials); -periodicSignal(end:-1:1,:)];
            else % If odd length
                periodicSignal = 2*pi*rand((obj.fftLength-1)/2,obj.trials);
                paddedPeriodicSignal = [zeros(1,obj.trials); periodicSignal;...
                    -periodicSignal(end:-1:1,:)];
            end
            stimFreqDomain = exp(i* paddedPeriodicSignal);
            timeDomainSignal = ifft(stimFreqDomain,'symmetric');
            switch obj.arithmetic
                case 'single'
                    timeDomainSignal = single(timeDomainSignal);
                    stimFreqDomain = fft(timeDomainSignal);
                    stimTimeDomain = [timeDomainSignal; timeDomainSignal];
                case 'double'
                    timeDomainSignal = double(timeDomainSignal);
                    stimTimeDomain = [timeDomainSignal; timeDomainSignal];
                case 'signed'
                    % Scale numbers before casting in fixed point
                    r = 2^(obj.iwl-obj.ifl-1); % The "-1" is needed for onesided range
                    timeDomainSignalmax = max(max(abs(timeDomainSignal)));
                    timeDomainSignal = r*timeDomainSignal./timeDomainSignalmax;
                    stimFreqDomain = fft(timeDomainSignal);
                    stimTimeDomain = [timeDomainSignal; timeDomainSignal];
                    stimTimeDomain = sfi(stimTimeDomain,obj.iwl,obj.ifl);
                case 'unsigned'
                    % Scale numbers before casting in fixed point
                    r = 2^(obj.iwl-obj.ifl-1); % The "-1" is needed for onesided range
                    timeDomainSignalmax = max(max(abs(timeDomainSignal)));
                    timeDomainSignal = r*timeDomainSignal./timeDomainSignalmax;
                    stimFreqDomain = fft(timeDomainSignal);
                    stimTimeDomain = [timeDomainSignal; timeDomainSignal];
                    stimTimeDomain = ufi(stimTimeDomain,obj.iwl,obj.ifl);
                otherwise
                    stimTimeDomain = [timeDomainSignal; timeDomainSignal];
            end
            
        end
        function filterOutFreqDomain = getFilterResponse(obj, stimTimeDomain, myFilter)
            filterOutTimeDomain = myFilter(stimTimeDomain);
            filterOutTimeDomain = double(filterOutTimeDomain);
            filterOutTimeDomain = filterOutTimeDomain(obj.fftLength+1:2*obj.fftLength,:);
            filterOutFreqDomain = fft(filterOutTimeDomain);
        end
        function NoisePowerSpectralDensity = theMagic(obj, filterOutFreqDomain, stimFreqDomain)
            sumFilterOutFreqDomainTrials = sum(filterOutFreqDomain.*conj(filterOutFreqDomain),2);
            sumFilterTransferFunction = sum(filterOutFreqDomain./stimFreqDomain,2);
            filterTransferFunction = sumFilterTransferFunction/obj.trials;
            NoisePowerSpectralDensity = ((sumFilterOutFreqDomainTrials/(obj.trials-1)) - max(max(abs(stimFreqDomain)))^2*real(filterTransferFunction.*conj(filterTransferFunction))*obj.trials/(obj.trials-1))/obj.fftLength;
            switch obj.arithmetic
                case 'signed'
                    % Later Implement Roundmode Fix makes d = 3;
                    d = 12;
                    R = 2^(obj.iwl-obj.ifl); % Input Quantizer Range
                    delta = R/(2^obj.iwl);
                    errvar = delta^2/d;
                    NoisePowerSpectralDensity = NoisePowerSpectralDensity - filterTransferFunction.*conj(filterTransferFunction)*errvar;
                case 'unsigned'
                case 'single'
                otherwise
            end
            NoisePowerSpectralDensity = abs(NoisePowerSpectralDensity); % Make sure it is nonnegative despite roundoff
        end
    end
end
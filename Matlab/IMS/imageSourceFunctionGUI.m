% imageSourceFunction.m

% %%Init Matlav
% clear all;
% close all;
function IMS(handles)
%%
%%variables
Lx = get(handles.RSX, 'Value'); %Length
Ly = get(handles.RSX, 'Value'); %Width
Lz = get(handles.RSX, 'Value'); %Height
a = get(handles.LLX, 'Value');  %Listener loc X
b = get(handles.LLY, 'Value'); %Listener loc Y
c = get(handles.LLZ, 'Value'); %Listener loc Z
p = get(handles.SLX, 'Value'); %Source loc X
q = get(handles.SLY, 'Value'); %Source loc Y
r = get(handles.SLZ, 'Value'); %Source loc Z
alphaxpos = get(handles.alphaxpos, 'Value'); %Reflection Coefficient
alphaxneg = get(handles.alphaxneg, 'Value'); %Reflection Coefficient
alphaypos = get(handles.alphaypos, 'Value'); %Reflection Coefficient
alphayneg = get(handles.alphayneg, 'Value'); %Reflection Coefficient
alphazpos = get(handles.alphazpos, 'Value'); %Reflection Coefficient
alphazneg = get(handles.alphazneg, 'Value'); %Reflection Coefficient
cs = 343; %Speed of sound
% [audio Fs] = audioread('it could be sweet.mp3'); %Sampling Frequency
N = get(handles.MRO, 'Value'); % Order of reflections
% lpfco = 8000; %low pass filter cut off
% filtord = 6; %filter order#
dE = 0.2; %distance between ears
% astart = 20;%start of audio to process (s)
% afinish = 40;%end of audio to process (s)
%%
%Calc sum stuff


%Calc difference between arival times for ears 
 theta = atan((p-a)/(q-b));

%Calc for the ear positions 
aleft = a + (dE/2) * cos(theta);
aright = a - (dE/2) * cos(theta);

bleft = b - (dE/2) * sin(theta);
bright = b + (dE/2) * sin(theta);

%calc sim time
Tf = (1/cs) * sqrt((N*Lx+p-a)^2+(N*Ly+q-b)^2+(N*Lz+r-c)^2);


%Initz IR Vector
irl = zeros(1,ceil(Tf * Fs));
irr = zeros(1,ceil(Tf * Fs));
%%
%Do Processing for left

for d = -N : N
    %Find Distance
    if (mod(d,2) == 1)
        % if odd
        A = (d+1)* Lx - p - aleft;
    else
        % if even
        A = d * Lx + p - aleft;
    end
    
    for e = -N : N
        %Find Distance
        if (mod(e,2) == 1)
            % if odd
            B = (e+1)* Ly - q - bleft;
        else
            % if even
            B = e* Ly + q - bleft;
        end
        
        for f = -N : N
            %Find Distance
            if (mod(f,2) == 1)
                % if odd
                C = (f+1) *Lz - r - c;
            else
                % if even
                C = f *Lz  + r - c;
            end
            %Calculate number of passes through each wall
            %LR
            [xposhits xneghits] = wallHits(d);
            %FB
            [yposhits yneghits] = wallHits(e);
            %UD
            [zposhits zneghits] = wallHits(f);
            
            %Calculate Distances
            L = sqrt((A^2)+(B^2)+(C^2));
            
            %Calculate reflection for each axis
            gA = (alphaxneg^(xneghits)) * (alphaxpos^(xposhits));
            
            gB = (alphayneg^(yneghits)) * (alphaypos^(yposhits));
            
            gC = (alphazneg^(zneghits)) * (alphazpos^(zposhits));
            
            %Convert total Distance to time
            t = (1/cs)*L;
            
            %calc total G
            g = (gA * gB * gC) / L;
            
            %Convert for IR sampling rate
            ind = ceil(t * Fs);
            
            %Write to IR
            irl(1, ind) =  g;

        end
    end
end

%%
%Do Processing for right

for d = -N : N
    %Find Distance
    if (mod(d,2) == 1)
        % if odd
        A = (d+1)* Lx - p - aright;
    else
        % if even
        A = d * Lx + p - aright;
    end
    
    for e = -N : N
        %Find Distance
        if (mod(e,2) == 1)
            % if odd
            B = (e+1)* Ly - q - bright;
        else
            % if even
            B = e* Ly + q - bright;
        end
        
        for f = -N : N
            %Find Distance
            if (mod(f,2) == 1)
                % if odd
                C = (f+1) *Lz - r - c;
            else
                % if even
                C = f *Lz  + r - c;
            end
            %Calculate number of passes through each wall
            %LR
            [xposhits xneghits] = wallHits(d);
            %FB
            [yposhits yneghits] = wallHits(e);
            %UD
            [zposhits zneghits] = wallHits(f);
            
            %Calculate Distances
            L = sqrt((A^2)+(B^2)+(C^2));
            
            %Calculate reflection for each axis
            gA = (alphaxneg^(xneghits)) * (alphaxpos^(xposhits));
            
            gB = (alphayneg^(yneghits)) * (alphaypos^(yposhits));
            
            gC = (alphazneg^(zneghits)) * (alphazpos^(zposhits));
            
            %Convert total Distance to time
            t = (1/cs)*L;
            
            %calc total G
            g = (gA * gB * gC) / L;
            
            %Convert for IR sampling rate
            ind = ceil(t * Fs);
            
            %Write to IR
            irr(1, ind) =  g;
        end
    end
end
set(handles.irl, 'Value', irl);
set(handles.irr, 'Value', irr);
% %combine tracks
% if length(irl) > length(irr)
% audioL = audioL(1:length(i));
% else
% audioR = audioR(1:length(audioL)); 
% end
% audioout(:,1) = audioL;
% audioout(:,2) = audioR;
% %Plot IR to make sure no silly
% plot(irl);
% hold on;
% plot(irr);
% hold off;
% %Calculate Filter 
% [B1, A1] = butter(filtord, lpfco/(Fs/2),'low');
% %Do Filtering for HF loss due to distance
% irl = filter(B1, A1, irl);
% irr = filter(B1, A1, irr);
% %Load Audio
% [audio afs] = audioread('it could be sweet.mp3');
% %Split audio
% audioL = audio((astart*Fs):(afinish*Fs),1);
% audioR = audio((astart*Fs):(afinish*Fs),2);
% %Do Conv
% audioL = conv(audioL,irl);
% audioR = conv(audioR,irr);
% %combine tracks
% if length(audioL) > length(audioR)
% audioL = audioL(1:length(audioR));
% else
% audioR = audioR(1:length(audioL)); 
% end
% audioout(:,1) = audioL;
% audioout(:,2) = audioR;
% %Normalise Conv
% audioout = audioout/(max(abs(audioout)));
% %Play and compare
% player2 = audioplayer(audio((astart*Fs):(afinish*Fs),:),Fs);
% player = audioplayer(audioout, Fs);
% playblocking(player2);
% play(player);
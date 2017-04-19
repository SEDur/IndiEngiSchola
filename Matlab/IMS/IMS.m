% imageSourceFunction.m

% %%Init Matlav
% clear all;
% close all;
function handles = IMS(handles)
%%
%%variables
Lx = get(handles.RSX, 'Value'); %Length
Ly = get(handles.RSX, 'Value'); %Width
Lz = get(handles.RSX, 'Value'); %Height
a = get(handles.LLX, 'Value');  %Listener loc X
b = get(handles.LLY, 'Value'); %Listener loc Y
c = get(handles.LLZ, 'Value'); %Listener loc Z
% Source 1 infro
p1 = get(handles.SLX1, 'Value'); %Source loc X
q1 = get(handles.SLY1, 'Value'); %Source loc Y
r1 = get(handles.SLZ1, 'Value'); %Source loc Z
% Source 2 infro
p2 = get(handles.SLX2, 'Value'); %Source loc X
q2 = get(handles.SLY2, 'Value'); %Source loc Y
r2 = get(handles.SLZ2, 'Value'); %Source loc Z

alphaxposs = handles.alphaxpos; %Reflection Coefficient
alphaxnegs = handles.alphaxneg; %Reflection Coefficient
alphayposs = handles.alphaypos; %Reflection Coefficient
alphaynegs = handles.alphayneg; %Reflection Coefficient
alphazposs = handles.alphazpos; %Reflection Coefficient
alphaznegs = handles.alphazneg; %Reflection Coefficient
cs = 343; %Speed of sound
Fs = handles.Fs;
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
 theta1 = atan((p1-a)/(q1-b));
 theta2 = atan((p2-a)/(q2-b));
%Calc for the ear positions 
aleft1 = a + (dE/2) * cos(theta1);
aright1 = a - (dE/2) * cos(theta1);

bleft1 = b - (dE/2) * sin(theta1);
bright1 = b + (dE/2) * sin(theta1);

aleft2 = a + (dE/2) * cos(theta2);
aright2 = a - (dE/2) * cos(theta2);

bleft2 = b - (dE/2) * sin(theta2);
bright2 = b + (dE/2) * sin(theta2);

%calc sim time
Tf1 = (1/cs) * sqrt((N*Lx+p1-a)^2+(N*Ly+q1-b)^2+(N*Lz+r1-c)^2);
Tf2 = (1/cs) * sqrt((N*Lx+p2-a)^2+(N*Ly+q2-b)^2+(N*Lz+r2-c)^2);

%Initz IR Vector
irl1 = zeros(1,ceil(Tf1 * Fs));
irr1 = zeros(1,ceil(Tf1 * Fs));
irl2 = zeros(1,ceil(Tf2 * Fs));
irr2 = zeros(1,ceil(Tf2 * Fs));
%%
%Do Processing for left

for d = -N : N
    %Find Distance
    if (mod(d,2) == 1)
        % if odd
        A1 = (d+1)* Lx - p1 - aleft1;
        A2 = (d+1)* Lx - p2 - aleft2;
    else
        % if even
        A1 = d * Lx + p1 - aleft1;
        A2 = d * Lx + p2 - aleft2;
    end
    
    for e = -N : N
        %Find Distance
        if (mod(e,2) == 1)
            % if odd
            B1 = (e+1)* Ly - q1 - bleft1;
            B2 = (e+1)* Ly - q2 - bleft2;
        else
            % if even
            B1 = e* Ly + q1 - bleft1;
            B2 = e* Ly + q2 - bleft2;
        end
        
        for f = -N : N
            %Find Distance
            if (mod(f,2) == 1)
                % if odd
                C1 = (f+1) *Lz - r1 - c;
                C2 = (f+1) *Lz - r2 - c;
            else
                % if even
                C1 = f *Lz  + r1 - c;
                C2 = f *Lz  + r2 - c;
            end
            %Calculate number of passes through each wall
            %LR
            [xposhits xneghits] = wallHits(d);
            %FB
            [yposhits yneghits] = wallHits(e);
            %UD
            [zposhits zneghits] = wallHits(f);
            
            %Calculate Distances
            L1 = sqrt((A1^2)+(B1^2)+(C1^2));
            L2 = sqrt((A2^2)+(B2^2)+(C2^2));
            %Calculate reflection for each axis
            gA = (alphaxnegs^(xneghits)) * (alphaxposs^(xposhits));
            
            gB = (alphaynegs^(yneghits)) * (alphayposs^(yposhits));
            
            gC = (alphaznegs^(zneghits)) * (alphazposs^(zposhits));
            
            %Convert total Distance to time
            t1 = (1/cs)*L1;
            t2 = (1/cs)*L2;
            %calc total G
            g1 = (gA * gB * gC) / L1;
            g2 = (gA * gB * gC) / L2;
            %Convert for IR sampling rate
            ind1 = ceil(t1 * Fs);
            ind2 = ceil(t2 * Fs);
            %Write to IR
            irl1(1, ind1) =  g1;
            irl2(1, ind2) =  g2;

        end
    end
end

%%
%Do Processing for right

for d = -N : N
    %Find Distance
    if (mod(d,2) == 1)
        % if odd
        A1 = (d+1)* Lx - p1 - aright1;
        A2 = (d+1)* Lx - p2 - aright2;
    else
        % if even
        A1 = d * Lx + p1 - aright1;
        A2 = d * Lx + p2 - aright2;
    end
    
    for e = -N : N
        %Find Distance
        if (mod(e,2) == 1)
            % if odd
            B1 = (e+1)* Ly - q1 - bright1;
            B2 = (e+1)* Ly - q2 - bright2;
        else
            % if even
            B1 = e* Ly + q1 - bright1;
            B2 = e* Ly + q2 - bright2;
        end
        
        for f = -N : N
            %Find Distance
            if (mod(f,2) == 1)
                % if odd
                C1 = (f+1) *Lz - r1 - c;
                C2 = (f+1) *Lz - r2 - c;
            else
                % if even
                C1 = f *Lz  + r1 - c;
                C2 = f *Lz  + r2 - c;
            end
            %Calculate number of passes through each wall
            %LR
            [xposhits xneghits] = wallHits(d);
            %FB
            [yposhits yneghits] = wallHits(e);
            %UD
            [zposhits zneghits] = wallHits(f);
            
            %Calculate Distances
            L1 = sqrt((A1^2)+(B1^2)+(C1^2));
            L2 = sqrt((A2^2)+(B2^2)+(C2^2));
            %Calculate reflection for each axis
            gA = (alphaxnegs^(xneghits)) * (alphaxposs^(xposhits));
            
            gB = (alphaynegs^(yneghits)) * (alphayposs^(yposhits));
            
            gC = (alphaznegs^(zneghits)) * (alphazposs^(zposhits));
            
            %Convert total Distance to time
            t1 = (1/cs)*L1;
            t2 = (1/cs)*L2;
            %calc total G
            g1 = (gA * gB * gC) / L1;
            g2 = (gA * gB * gC) / L2;
            %Convert for IR sampling rate
            ind1 = ceil(t1 * Fs);
            ind2 = ceil(t2 * Fs);
            %Write to IR
            irr1(1, ind1) =  g1;
            irr2(1, ind2) =  g2;
        end
    end
end
% set(handles.irl, 'Value', irl);
% set(handles.irr, 'Value', irr);

if length(irl1) > length(irl2)
irl1 = irl1(1:length(irl2));
else
irl2 = irl2(1:length(irl1)); 
end

if length(irr1) > length(irr2)
irr1 = irr1(1:length(irr2));
else
irr2 = irr2(1:length(irr1)); 
end


handles.irl = (irl1 + irl2);
handles.irr = (irr1 + irr2);
listen = 1;
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
end
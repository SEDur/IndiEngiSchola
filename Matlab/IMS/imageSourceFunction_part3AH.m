% imageSourceFunction
% Adapted from Adam Hills method
% Simon Durbridge
% Feb 18th 2016

clear all;
close all;


Fs = 48000;
N = 10;
Lx = 12;
Ly = 21;
Lz = 3;
a = 15;
b = 21;
c = 3;
p = 30;
q = 22;
r = 12;
alphaXpos = 0.9;
alphaXneg = 0.1;
alphaYpos = 0.9;
alphaYneg = 0.7;
alphaZpos = 0.1;
alphaZneg = 0.4;
cs = 343;
dE = 0.2;

% Calculate the two listener locations. They will be equal distances from
% the sound source.
theta = atan((p-a)/(q-b));

% calculate the maximum distance that tF can take on
tF = (1/cs)*sqrt((N*Lx + p - a)^2 + (N*Ly + q - b)^2 + (N*Lz + r - c)^2);

% Calculate the location of each of the user's "ears"
a_R = a - (dE/2)*cos(theta);
a_L = a + (dE/2)*cos(theta);
b_R = a + (dE/2)*sin(theta);
b_L = a - (dE/2)*sin(theta);

% initialize an impulse response vector with appropriate length
ir_R = zeros(1, ceil(tF*Fs));
ir_L = zeros(1, ceil(tF*Fs));

% *** RIGHT EAR ***
% Nested for loops to calculate all distances, amplitudes and then place
% the proper value in the iir vector (only deal with one at a time, no
% values are stored)
for d = -N : N                              % loop for x dimension
    % determine if d is even or odd
    if mod(d, 2) == 1
        A = (d + 1)*Lx - p - a_R;
    else
        A = d*Lx + p - a_R;
    end
    for e = -N : N                          % loop for y dimension
        % determine if e is even or odd
        if mod(e, 2) == 1
            B = (e + 1)*Ly - q - b_R;
        else
            B = e*Ly + q - b_R;
        end
        for f = -N : N                      % loop for z dimension
            % determine if f is even or odd
            if mod(f, 2) == 1
                C = (f + 1)*Lz - r - c;
            else
                C = f*Lz + r - c;
            end
            %**************************************************************
            % Find the number of times each wall is encountered
            %   Note that front + back hits shoulf equal d, left + right
            %   hits sould equal e, and floor + ceiling hits should equal f
            %**************************************************************
            [XNeg Xpos] = wallHits(d);
            [Yneg Ypos] = wallHits(e);
            [Zneg Zpos] = wallHits(f); 
            
            % cacluate the reflection from the vectors coming through the
            % front/back walls
            gA = (alphaXneg^(XNeg))*(alphaXpos^(Xpos));
            % cacluate the reflection from the vectors comming through the
            % left/right walls
            gB = (alphaYneg^(Yneg))*(alphaYpos^(Ypos));
            % cacluate the reflection from the vectors comming through the
            % floor and ceiling
            gC = (alphaZneg^(Zneg))*(alphaZpos^(Zpos));
            % calculate the length the impulse has traveled
            L = sqrt(A^2 + B^2 + C^2);
            % calculate the time the impulse arrives at the listener
            t = (1/343)*L;
            % calculate the amplitude the impulse has at the listener 
            g = (gA*gB*gC)/L;
            
            % add the impulse to the ir vector
            ind = ceil(t*Fs);
            ir_R(1, ind) = g;
        end
    end
end

% *** LEFT EAR ***
% Nested for loops to calculate all distances, amplitudes and then place
% the proper value in the iir vector (only deal with one at a time, no
% values are stored)
for d = -N : N                              % loop for x dimension
    % determine if d is even or odd
    if mod(d, 2) == 1
        A = (d + 1)*Lx - p - a_L;
    else
        A = d*Lx + p - a_L;
    end
    for e = -N : N                          % loop for y dimension
        % determine if e is even or odd
        if mod(e, 2) == 1
            B = (e + 1)*Ly - q - b_L;
        else
            B = e*Ly + q - b_L;
        end
        for f = -N : N                      % loop for z dimension
            % determine if f is even or odd
            if mod(f, 2) == 1
                C = (f + 1)*Lz - r - c;
            else
                C = f*Lz + r - c;
            end
            %**************************************************************
            % Find the number of times each wall is encountered
            %   Note that front + back hits shoulf equal d, left + right
            %   hits sould equal e, and floor + ceiling hits should equal f
            %**************************************************************
            [XNeg Xpos] = wallHits(d);
            [Yneg Ypos] = wallHits(e);
            [Zneg Zpos] = wallHits(f); 
            
            % cacluate the reflection from the vectors coming through the
            % front/back walls
            gA = (alphaXneg^(XNeg))*(alphaXpos^(Xpos));
            % cacluate the reflection from the vectors comming through the
            % left/right walls
            gB = (alphaYneg^(Yneg))*(alphaYpos^(Ypos));
            % cacluate the reflection from the vectors comming through the
            % floor and ceiling
            gC = (alphaZneg^(Zneg))*(alphaZpos^(Zpos));
            % calculate the length the impulse has traveled
            L = sqrt(A^2 + B^2 + C^2);
            % calculate the time the impulse arrives at the listener
            t = (1/343)*L;
            % calculate the amplitude the impulse has at the listener 
            g = (gA*gB*gC)/L;
            
            % add the impulse to the ir vector
            ind = ceil(t*Fs);
            ir_L(1, ind) = g;
        end
    end
end


% Apply a low pass filter to ir to simulate loss of high
% frequencies over time as the room absorbs the impulses (cuttoff: 8000Hz)
[B1, A1] = butter(6, 8000/(Fs/2), 'low');
ir_R = filter(B1, A1, ir_R);
ir_L = filter(B1, A1, ir_L);

% load in .wav file
% [dry, Fs] = audioread('dry.wav');
[predry Fs] = audioread('it could be sweet.mp3');

dry = predry( ((Fs*10)+1) : ((Fs*40)+1),:);

% normalize the dry signal
dry = dry/max(abs(dry));

% convolve dry signal with IR
wet_R = conv(dry, ir_R);
wet_L = conv(dry, ir_L);

% make stereo signal
if length(wet_L) > length(wet_R)
    wet_L = wet_L(1 : length(wet_R));
else
    wet_R = wet_R(1 : length(wet_L));
end
wet = [wet_L; wet_R];

% normalize wet signal
wet = wet/max(max(abs(wet)));

% playback dry sound
player1 = audioplayer(dry, Fs);
playblocking(player1);

% playback dry sound
player2 = audioplayer(wet, Fs);
play(player2);
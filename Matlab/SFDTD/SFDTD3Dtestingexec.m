%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SFDTD3Dtestingeec.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% A script that was written to test the execution time of the SFDTD method
% for domains of five different sizes.
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Initz Matlab
%%Initz Variables
%Units

%Distance
meters      = 1;
centimeters = 1e-2 * meters;
millimeters = 1e-3 * meters;
inches      = 2.54 * centimeters;
feet        = 12 * inches;
%Time
seconds     = 1; 
hertz       = 1/seconds;
kilohertz   = 1e3 * hertz;
megahertz   = 1e6 * hertz;
gigahertz   = 1e9 * hertz;

% constants
c     = 343 * meters / seconds; %Speed of sound m/s
rho    = 1.21; %Density of air kg/m^3
p0 = 2*10^-5;
cstab = 2/(pi*sqrt(3));
%% Initz Winodws
figure(1);
set(1, 'WindowStyle', 'Docked')
figure(2);
set(2, 'WindowStyle', 'Docked')

%% Hard Coded Variables

%Maximum calculation frequency
fmax = 500 * hertz;

% %grid size
gx = (c / (6*fmax));
gy = (c / (6*fmax));
gz = (c / (6*fmax));
dt = ((1/c)*gx)/2;

%Boundary Absorption Coefs (0 to 1)
alphaL = 0.45;
alphaR = 0.45;
alphaF = 0.45;
alphaB = 0.45;
alphaT = 0.45;
alphaG = 0.45;

%number of sources
snum = 1;

%Time of sim
T = 1.0 ;
% generate the source(s) & determine number of time steps needed

%% Gauss Source
% tnum = ceil(T/dt);
% fc = 0.05;     % Cutoff frequency (normalised 0.5=nyquist)
% n0 = 30;        % Initial delay (samples)
% sigma=sqrt(2*log(2))/(2*pi*(fc/dt));
% n=0:tnum;
% source1=exp(-dt^2*(n-n0).^2/(2*sigma^2)).*(10^-12*10^(80/20));
% for n = 37 : length(source1)
%     if(source1(n) < 0)
%        source1(n) = 0; 
%     end
% end
%% Toneburst Source
% tone = dsp.SineWave('Amplitude',((2*10^-5)*10^(100/20)),...
%     'Frequency', 1000,...
%     'SampleRate', 1/dt,...
%     'SamplesPerFrame',0.01/dt);
% w1 = window(@gausswin,0.01/dt,2.5); 
% toneBurst = tone() .* w1;
% source1 = zeros(T/dt,1);
% source1(10:69) = toneBurst;
% source1(110:169) = toneBurst;
% source1(410:469) = toneBurst;

%% MLS Source
source1 = GenerateMLSSequence(5,9,0).*((2*10^-5)*10^(100/20));
T = length(source1)*dt;
w1 = window(@gausswin,length(source1),2.5); 
source1 = source1 .* w1;

widths = [5 10 20 40 60];
for cntr = 1 : 5
%define grid width
lx = widths(cntr);
ly = lx;
lz = lx;

%number of cells
xcells = ceil(lx/gx);
ycells = ceil(ly/gy);
zcells = ceil(lz/gz);

%source locations
sourcelocations = [ceil(1/gy) ceil(1/gx) ceil(1/gz);...
                    ceil(xcells/2) ceil(ycells/2) ceil(zcells/2)];
%% Pre cast pressure and velocity matrices
% initialize the velocity and pressure matrices (matrices are set up in a
% y by x fashion to properly display the 2D space (y = rows, x = columns))
p = ones(ycells - 1, xcells - 1, zcells - 1)*(2*10^-5);
ux = zeros(ycells - 1, xcells, zcells - 1);
uy = zeros(ycells, xcells - 1, zcells - 1);
uz = zeros(ycells - 1, xcells - 1, zcells);

%% Calculate coefficients
% set up the multiplication constants for the update equations
uCx = dt/(gx*rho);
uCy = dt/(gy*rho);
uCz = dt/(gz*rho);
pCx = c^2*rho*dt/gx;
pCy = c^2*rho*dt/gy;
pCz = c^2*rho*dt/gz;

% set the characteristic impedances of the walls
ZR = rho*c*(1 + sqrt(1 - alphaR))/(1 - sqrt(1 - alphaR));
ZL = rho*c*(1 + sqrt(1 - alphaL))/(1 - sqrt(1 - alphaL));
ZF = rho*c*(1 + sqrt(1 - alphaF))/(1 - sqrt(1 - alphaF));
ZB = rho*c*(1 + sqrt(1 - alphaB))/(1 - sqrt(1 - alphaB));
ZT = rho*c*(1 + sqrt(1 - alphaT))/(1 - sqrt(1 - alphaT));
ZG = rho*c*(1 + sqrt(1 - alphaG))/(1 - sqrt(1 - alphaG));
% calulcate the coefficients used for the boundary conditions
Rx = rho*gx/dt;
Ry = rho*gy/dt;
Rz = rho*gz/dt;

% plot vectors
linex = linspace(0, lx - gx, xcells-1);
liney = linspace(0, ly - gy, ycells-1);
linez = linspace(0, lz - gz, zcells-1);

%Set zsliceloc
[xvec, yvec, zvec] = meshgrid(0 : gx : gx * (xcells-2),...
    0 : gy : gy * (ycells-2), 0 : gz : gz * (zcells-2));

%% Run for time
for n = 1:T/dt    
    disp(T-(n*dt))
    tic();
    [idx] = SPARSEfun3DC(p, 10, p0);
    [p, ux, uy, uz] = SFDTD3DfunC(p, pCx, pCy, pCz, ux, uy, uz, uCx, uCy, uCz, Rx, Ry, Rz, ZL, ZR, ZT, ZB, ZF, ZB, idx);
    reciever(n,1) = p(ceil(xcells/4), ceil(ycells/4),ceil(zcells/2));
    reciever(n,2) = p(ceil(xcells/2)+ceil(xcells/4), ceil(ycells/4),ceil(zcells/2));    
    reciever(n,3) = p(ceil(xcells/2), ceil(ycells/2)+ceil(xcells/4),ceil(zcells/2));
    reciever(n,4) = p(ceil(xcells/2+ceil(xcells/4)), ceil(ycells/2)+ceil(xcells/4),ceil(zcells/2));
    reciever(n,5) = p(ceil(xcells/2), ceil(ycells/2),ceil(zcells/2));
    srcnorm(n) = p(sourcelocations(1,1),sourcelocations(1,2),sourcelocations(1,3));
    % Input source
    p(sourcelocations(1),sourcelocations(2),sourcelocations(3)) = p(sourcelocations(1),sourcelocations(2),sourcelocations(3)) - source1(n);
    exectime(n) = toc();
end

%% Some really postprocessing
norec = reciever ./ max(abs(reciever));
% recanal = AnalyseMLSSequence(reciever',0,5,9,0,0);
[lpsd, lf] = pwelch(norec,hann(1000),[],1000,1/dt);
srcnrm = srcnorm ./ max(abs(srcnorm));
[spsd, sf] = pwelch(srcnrm,hann(1000),[],1000,1/dt);
filename = strcat('xwidth',num2str(lx),'.mat');
save(filename,'exectime', 'norec', 'lpsd', 'lf', 'srcnrm', 'spsd', 'sf');
end
%% Postprocessing

for i = 1 : size(reciever,2)
    lag(i) = getlag(reciever(:,i)',srcnorm);
    recieverCirc(:,i) = circshift(reciever(:,i)',lag(i));
end

% Hd = postprocessingDCfilter;
norec = recieverCirc ./ max(abs(recieverCirc));
% recanal = AnalyseMLSSequence(reciever',0,3,11,0,0);
% norec = Hd(norec);
% [lpsd, lf] = pwelch(norec,hann(5000),[],5000,fs);
[lpsd, lf] = pwelch(norec,hann(5000),[],5000,1/dt);
% clear('Hd');
% Hd = postprocessingDCfilter;
srcnrm = srcnorm ./ max(abs(srcnorm));
% srcnrm = Hd(srcnrm);
% [spsd, sf] = pwelch(srcnrm,hann(5000),[],5000,fs);
[spsd, sf] = pwelch(srcnrm,hann(5000),[],5000,1/dt);

%% Display the results
subplot(3,1,1);
plot(0:dt:((length(reciever)-1)*dt),source1(1:length(reciever)))
hold on;
plot(0:dt:((length(reciever)-1)*dt),reciever)
hold off;
axis('tight')
legend('source','rec top left','rec top right','rec bottom left',...
    'rec bottom right','rec centre');
title('Raw Input And Output');
xlim([0 0.2]);
xlabel('Time (s)');
ylabel('Amplitude (Pa)');
subplot(3,1,2);
plot(0:dt:((length(norec)-1)*dt),norec,'--','linewidth',2.0)
hold on;
plot(0:dt:((length(srcnrm)-1)*dt),srcnrm,'-')
hold off;
axis('tight')
legend('rec top left','rec top right','rec bottom left',...
    'rec bottom right','rec centre','source');
title('Normalised Input And Output');
xlim([0 0.2]);
xlabel('Time (s)');
ylabel('Amplitude (Pa)');
subplot(3,1,3);
plot(lf, db(lpsd),'--','Linewidth',2.0);
hold on;
plot(sf, db(spsd));
hold off;
legend('rec top left','rec top right','rec bottom left',...
    'rec bottom right','rec centre','source');
grid('on');
title('Power Spectral Density of Input and Output');
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');

% 3D implementation
% S Durbridge 
% 2017

%%Initz Matlab
% clear all;
% close all;
%Units
figure(1);
set(1, 'WindowStyle', 'Docked');
cstab = 2/(pi*sqrt(3));
meters = 1;
hertz = 1;
c = 343;
p0 = 2*10^-5;
rho = 1.21;

%%
%%Hard Code Variables
%Maximum calculation frequency
fmax = 5000 * hertz;
% dt = 1/ (c*sqrt((1/(gy^2))+(1/(gx^2))+(1/(gz^2))));
% dt = 1/(2*fmax);

%grid size
% gx = c * dt / cstab;
% gy = c * dt / cstab;
% gz = c * dt / cstab;
gx = (c / (6*fmax));
gy = (c / (6*fmax));
gz = (c / (6*fmax));
dt = ((1/c)*gx)/2;
%Dims
%Dim Size (m)
lx = 5*meters;
ly = 4*meters;
lz = 3*meters;

xcells = ceil(lx/gx);
ycells = ceil(ly/gy);
zcells = ceil(lz/gz);

%number of sources
snum = 1;
%source locations
% sourcelocations = [ceil((ly/gy)/2) ceil((lx/gx)/2) ceil((lz/gz)/2);...
%                     ceil((ly/gy)/2)+1 ceil((lx/gx)/2)+1 ceil((lz/gz)/2)+1];
sourcelocations = [ceil((1/gy)) ceil(1/gx) ceil(1/gz);...
                    ceil((1/gy)) ceil(1/gx) ceil(1/gz)];

%recieves position
recieverleftloc = [ceil(ycells/2) ceil(xcells/2) ceil(zcells/2)];
%Time of sim
T = 1.0;

% generate the source(s) & determine number of time steps needed
% tnum = ceil(T/dt);
% fc = 0.1;     % Cutoff frequency (normalised 0.5=nyquist)
% n0 = 30;        % Initial delay (samples)
% sigma=sqrt(2*log(2))/(2*pi*(fc/dt));
% n=0:tnum;
% source1=exp(-dt^2*(n-n0).^2/(2*sigma^2));
% source1= (source1 ./ max(source1)).*((2*10^-5)*10^(100/20));
% w1 = window(@hamming,0.4/(dt)); 
% fir = dsp.FIRFilter;
% fir.Numerator = w1';
% chirp = dsp.Chirp(...
%     'SweepDirection', 'Unidirectional', ...
%     'TargetFrequency', ceil(fmax/2), ...
%     'InitialFrequency', 10,...
%     'TargetTime', 0.4, ...
%     'SweepTime', 0.4, ...
%     'SamplesPerFrame', 0.4/dt, ...
%     'SampleRate', 1/dt);
% % plot(chirp());
% source1 = w1.*chirp();
% source1 = [zeros(10,1); source1];
% source1 = [source1; zeros((T/dt) - length(source1),1)].*((2*10^-5)*10^(100/20));

source1 = GenerateMLSSequence(3,11,0).*((2*10^-5)*10^(100/20));
T = length(source1)*dt;
w1 = window(@gausswin,length(source1),2.5); 
source1 = source1 .* w1;
% initialize the velocity and pressure matrices (matrices are set up in a
% y by x fashion to properly display the 2D space (y = rows, x = columns))
% p = ones(ycells - 1, xcells - 1, zcells - 1) .* 10^-12*10^(40/20);
p = zeros(ycells - 1, xcells - 1, zcells - 1);
ux = zeros(ycells - 1, xcells, zcells - 1);
uy = zeros(ycells, xcells - 1, zcells - 1);
uz = zeros(ycells - 1, xcells - 1, zcells);

% set up the multiplication constants for the update equations
uCx = dt/(gx*rho);
uCy = dt/(gy*rho);
uCz = dt/(gz*rho);
pCx = c^2*rho*dt/gx;
pCy = c^2*rho*dt/gy;
pCz = c^2*rho*dt/gz;

% set the wall reflection coefficients
% if alphaX = 0, then slightly adjust to avoid infinite characteristic
% impedance.

%Boundary Absorption Coefs (0 to 1)
alphaL = 0.45;
alphaR = 0.45;
alphaF = 0.45;
alphaB = 0.45;
alphaT = 0.45;
alphaG = 0.45;

if alphaR == 0
   alphaR = 1e-016; 
end
if alphaL == 0
   alphaL = 1e-016; 
end
if alphaF == 0
   alphaF = 1e-016; 
end
if alphaB == 0
   alphaB = 1e-016; 
end
if alphaT == 0
   alphaT = 1e-016; 
end
if alphaG == 0
   alphaG = 1e-016; 
end
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

%Initialize recording vectors
% leftear = zeros(1,tnum);
% rightear = zeros(1,tnum);
% meanpstore = zeros(1,tnum);

% xvec = 0 : gx : gx * (xcells-2);
% yvec = 0 : gy : gy * (ycells-2);
% zvec = 0 : gz : gz * (zcells-2);

[xvec, yvec, zvec] = meshgrid(0 : gx : gx * (xcells-2),...
    0 : gy : gy * (ycells-2), 0 : gz : gz * (zcells-2));
%%
%Set zsliceloc
% zslice = (s1loc(3));
% loop to update the velocities and pressures over the time steps, n
n = 1;
% while or((mean(mean(mean(abs(real(p(:,:,:)))))) > (p0 * 10^(60/10))),(n < (1600)))
while n*dt < T
    n = n + 1;
    T-(n*dt)
    tic();
    [p, ux, uy, uz] = FDTD3Dfun(p, pCx, pCy, pCz, ux, uy, uz, uCx,...
        uCy, uCz, Rx, Ry, Rz, ZL, ZR, ZF, ZB, ZT, ZG);
    p = FDTD3Dsources(p,sourcelocations ,source1(n) , 'soft');
    reciever(n) = p(recieverleftloc(1),recieverleftloc(2),recieverleftloc(3));
    srcnorm(n) = p(sourcelocations(1,1),sourcelocations(1,2),sourcelocations(1,3));
    exectime(n) = toc();
    FDTD3Dplotdomain(p, xcells, ycells, zcells, n, dt, p0); 
end
% % figure(2);
% subplot(3,1,1);
% plot(0:dt:((length(reciever)-1)*dt),reciever)
% axis('tight')
% subplot(3,1,2);
% plot(0:dt:((length(reciever)-1)*dt),source1)
% % xlim([0 0.01])
% subplot(3,1,3);
% plot(0:dt:((length(reciever)-1)*dt),exectime)
% axis('tight')

%% Some really minor postprocessing
% Hd = postprocessingDCfilter;
norec = reciever ./ max(abs(reciever));
% recanal = AnalyseMLSSequence(reciever',0,3,11,0,0);
% norec = Hd(norec);
% [lpsd, lf] = pwelch(norec,hann(5000),[],5000,fs);
[lpsd, lf] = pwelch(norec,hann(2000),[],2000,1/dt);
% clear('Hd');
% Hd = postprocessingDCfilter;
srcnrm = srcnorm ./ max(abs(srcnorm));
% srcnrm = Hd(srcnrm);
% [spsd, sf] = pwelch(srcnrm,hann(5000),[],5000,fs);
[spsd, sf] = pwelch(srcnrm,hann(2000),[],2000,1/dt);
%% Display the results
subplot(4,1,1);
plot(0:dt:((length(reciever)-1)*dt),reciever)
hold on;
plot(0:dt:((length(reciever)-1)*dt),source1(1:length(reciever)))
hold off;
axis('tight')
legend('reciever','source');
title('Raw Input And Output');
subplot(4,1,2);
plot(0:dt:((length(norec)-1)*dt),norec,'--','linewidth',2.0)
hold on;
plot(0:dt:((length(srcnrm)-1)*dt),srcnrm)
hold off;
axis('tight')
legend('reciever','source');
title('Normalised Input And Output');
subplot(4,1,3);
plot(lf, db(lpsd),'--','Linewidth',2.0);
hold on;
plot(sf, db(spsd));
hold off;
legend('reciever','source');
grid('on');
title('Power Spectral Density of Input and Output');
subplot(4,1,4);
plot(0:dt:((length(reciever)-1)*dt),roundtime)
axis('tight')
ttlstr = sprintf('Computation Time Per Cycle, Total Time: %i',sum(roundtime));
title(ttlstr);


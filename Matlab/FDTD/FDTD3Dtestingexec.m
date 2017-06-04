%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FDTD3Dtestingexec.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% A script that will test the execution speed of a 3D FDTD simulation with
% 5 different doman sizes.
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%Initz Matlab
% clear all;
% close all;
addpath(genpath('../mls/mls/'));
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
fmax = 500 * hertz;

%grid size
gx = (c / (6*fmax));
gy = (c / (6*fmax));
gz = (c / (6*fmax));
dt = ((1/c)*gx)/2;

%create source term
source1 = GenerateMLSSequence(2,7,0).*((2*10^-5)*10^(100/20));
w1 = window(@hamming,length(source1)); 
source1 = source1 .* w1;
T = length(source1)*dt;

% set up the multiplication constants for the update equations
uCx = dt/(gx*rho);
uCy = dt/(gy*rho);
uCz = dt/(gz*rho);
pCx = c^2*rho*dt/gx;
pCy = c^2*rho*dt/gy;
pCz = c^2*rho*dt/gz;

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

widths = [5 10 20 40 60];
for cntr = 1 : 5
%define grid width
lx = widths(cntr);
ly = lx;
lz = ly;

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

%preallocate domains
p = zeros(ycells - 1, xcells - 1, zcells - 1);
ux = zeros(ycells - 1, xcells, zcells - 1);
uy = zeros(ycells, xcells - 1, zcells - 1);
uz = zeros(ycells - 1, xcells - 1, zcells);

% plot vectors
linex = linspace(0, lx - gx, xcells-1);
liney = linspace(0, ly - gy, ycells-1);
linez = linspace(0, lz - gz, zcells-1);

[xvec, yvec, zvec] = meshgrid(0 : gx : gx * (xcells-2),...
    0 : gy : gy * (ycells-2), 0 : gz : gz * (zcells-2));

%% loop to update the velocities and pressures over the time steps, n
n = 1;
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
%     FDTD3Dplotdomain(p, xcells, ycells, zcells, n, dt, p0); 
end
%% Postprocessing
% Hd = postprocessingDCfilter;
norec = reciever ./ max(abs(reciever));
recanal = AnalyseMLSSequence(reciever',0,2,5,0,0);
% norec = Hd(norec);
% [lpsd, lf] = pwelch(norec,hann(5000),[],5000,fs);
[lpsd, lf] = pwelch(norec,hann(ceil(length(norec)/2)),[],ceil(length(norec)/2),1/dt);
% clear('Hd');
% Hd = postprocessingDCfilter;
srcnrm = srcnorm ./ max(abs(srcnorm));
% srcnrm = Hd(srcnrm);
% [spsd, sf] = pwelch(srcnrm,hann(5000),[],5000,fs);
[spsd, sf] = pwelch(srcnrm,hann(ceil(length(norec)/2)),[],ceil(length(norec)/2),1/dt);
filename = strcat('xwidth',num2str(lx),'.mat');
save(filename,'exectime', 'norec', 'recanal', 'lpsd', 'lf', 'srcnrm', 'spsd', 'sf');
end
%% Display the results
subplot(5,1,1);
plot(0:dt:((length(reciever)-1)*dt),reciever)
hold on;
plot(0:dt:((length(reciever)-1)*dt),source1(1:length(reciever)))
hold off;
axis('tight')
legend('reciever','source');
title('raw input and output');
subplot(5,1,2);
plot(0:dt:((length(norec)-1)*dt),norec,'--','linewidth',2.0)
hold on;
plot(0:dt:((length(srcnrm)-1)*dt),srcnrm)
hold off;
axis('tight')
legend('reciever','source');
title('normalised input and output');
subplot(5,1,3);
plot(lf, db(lpsd));
hold on;
plot(sf, db(spsd));
hold off;
legend('reciever','source');
grid('on');
title('power spectral density of input and output');
subplot(5,1,4);
plot(0:dt:((length(reciever)-1)*dt),roundtime)
axis('tight')
ttlstr = sprintf('computation time per cycle, total time is %d',sum(roundtime));
title(ttlstr);
subplot(5,1,5);
plot(0:dt:((length(recanal)-1)*dt),recanal);
title('MLS Analysed');
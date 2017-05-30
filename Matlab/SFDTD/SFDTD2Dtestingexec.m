%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SFDTD2Dtestingexec.m
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

%% Initz Matlab
clf;
addpath(genpath('../mls/mls/'));
figure(3)
set(3, 'windowstyle','docked','color', 'w');

%% Initz Variables

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
cstab = 2/(pi*sqrt(2));
%%Hard Code Variables

%Maximum calculation frequency
fmax = 1000 * hertz;

%calculate space and time step
dx = (c / (6*fmax));
dy = (c / (6*fmax));
dt = ((1/c)*dx)/2;

%Boundary Absorption Coefs (0 to 1)
alphaL = 0.45;
alphaR = 0.45;
alphaF = 0.45;
alphaB = 0.45;

%Calculate coefficients for pressure and velocity
uCx = dt/(dx*rho);
uCy = dt/(dy*rho);
pCx = c^2*rho*dt/dx;
pCy = c^2*rho*dt/dy;

pidxRow = [];
pidxCol = [];
uxidx = [];
uyidx = [];

% set the wall reflection coefficients
% if alphaX = 0, then slightly adjust to avoid infinite characteristic
% impedance.
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
% set the characteristic impedances of the walls
ZR = rho*c*(1 + sqrt(1 - alphaR))/(1 - sqrt(1 - alphaR));
ZL = rho*c*(1 + sqrt(1 - alphaL))/(1 - sqrt(1 - alphaL));
ZT = rho*c*(1 + sqrt(1 - alphaF))/(1 - sqrt(1 - alphaF));
ZB = rho*c*(1 + sqrt(1 - alphaB))/(1 - sqrt(1 - alphaB));

% calulcate the coefficients used for the boundary conditions
Rx = rho*dx/dt;
Ry = rho*dy/dt;

source1 = GenerateMLSSequence(5,9,0).*((2*10^-5)*10^(100/20));
w1 = window(@hamming,length(source1)); 
source1 = source1 .* w1;
T = length(source1)*dt;

widths = [5 10 20 40 60];
for cntr = 1 : 5
%define grid width
lx = widths(cntr);
ly = lx;

%Dims
xcells = ceil(lx/dx);
ycells = ceil(lx/dy);

%number of sources
snum = 1;

%source locations
sourcelocations = [ceil((ycells/2)) ceil(xcells/2)];

%recieves position
recieverleftloc = [ceil(ycells/2) ceil(xcells/4)];

% initialize the velocity and pressure matrices (matrices are set up in a
% y by x fashion to properly display the 2D space (y = rows, x = columns))
p = zeros(ycells - 1, xcells - 1);
ux = zeros(ycells - 1, xcells);
uy = zeros(ycells, xcells - 1);
idx3 = uy;
% set up the multiplication constants for the update equations

% plot vectors
linex = linspace(0, lx - dx, xcells-1);
liney = linspace(0, ly - dy, ycells-1);

[xvec, yvec] = meshgrid(0 : dx : dx * (xcells-2),...
    0 : dy : dy * (ycells-2));

% loop to update the velocities and pressures over the time steps, n
n = 0;
for n = 1:T/dt    
n = n + 1;
    tic;
    [idx] = SPARSEfun2DC(p, 40, p0);
    [p, ux, uy] = SFDTD2Dfun(p, pCx, pCy, ux, uy, uCx, uCy, Rx, Ry, ZL,...
        ZR, ZT, ZB, idx);
    extime(n) = toc;
    p(sourcelocations(1),sourcelocations(2)) = p(sourcelocations(1),sourcelocations(2)) - source1(n);
    reciever(n) = p(recieverleftloc(1),recieverleftloc(2));
    srcnorm(n) = p(sourcelocations(1,1),sourcelocations(1,2));
    exectime(n) = toc();
end

%% Some really postprocessing
norec = reciever ./ max(abs(reciever));
recanal = AnalyseMLSSequence(reciever',0,5,9,0,0);
[lpsd, lf] = pwelch(norec,hann(1000),[],1000,1/dt);
srcnrm = srcnorm ./ max(abs(srcnorm));
[spsd, sf] = pwelch(srcnrm,hann(1000),[],1000,1/dt);
filename = strcat('xwidth',num2str(lx),'.mat');
save(filename,'exectime', 'norec', 'recanal', 'lpsd', 'lf', 'srcnrm', 'spsd', 'sf');
end
%% Display the results
subplot(4,1,1);
plot(0:dt:((length(reciever)-1)*dt),reciever)
hold on;
plot(0:dt:((length(reciever)-1)*dt),source1(1:length(reciever)))
hold off;
axis('tight')
legend('reciever','source');
title('raw input and output');
subplot(4,1,2);
plot(0:dt:((length(norec)-1)*dt),norec,'--','linewidth',2.0)
hold on;
plot(0:dt:((length(srcnrm)-1)*dt),srcnrm)
hold off;
axis('tight')
legend('reciever','source');
title('normalised input and output');
subplot(4,1,3);
plot(lf, db(lpsd),'--','Linewidth',2.0);
hold on;
plot(sf, db(spsd));
hold off;
legend('reciever','source');
grid('on');
title('power spectral density of input and output');
subplot(4,1,4);
plot(0:dt:((length(reciever)-1)*dt),exectime)
axis('tight')
ttlstr = sprintf('computation time per cycle, total time is %d',sum(exectime));
title(ttlstr);


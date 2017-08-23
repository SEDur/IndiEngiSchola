%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FDTD2Dtesting.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% A script that will execute an FDTD simulation in 2D.
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initz Matlab
clear all;
% close all;
addpath(genpath('../mls/mls/'));
figure(1)
set(1, 'windowstyle','docked','color', 'w');

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
% cstab = sqrt(1/2);
cstab = 2/(pi*sqrt(2));
% cstab = 1;

%%
%%Hard Code Variables
%Maximum calculation frequency
fmax = 5000 * hertz;

%grid size
gx = c * (1/fmax) / cstab;
gy = c * (1/fmax) / cstab;

%Dims
%Dim Size (m)
lx = 5*meters;
ly = 4*meters;
xcells = ceil(lx/gx);
ycells = ceil(lx/gy);

%Boundary Absorption Coefs (0 to 1)
alphaL = 0.45;
alphaR = 0.45;
alphaF = 0.45;
alphaB = 0.45;

%number of sources
snum = 2;
%source locations
s1loc = [ceil((lx/gx)/2) ceil((lx/gx)/2)];
s2loc = [ceil((ly/gx)/4) ceil((ly/gy)/2)];%source frequency

%recieves position
recieverleftloc = [ceil(ycells/2) ceil(xcells/2)];

%Time of sim
dt = 1/ (c*sqrt((1/(gy^2))+(1/(gx^2))));
T = 1;

% generate the source(s) & determine number of time steps needed
tnum = ceil(T/dt);
source1 = zeros(1,tnum);
source2 = zeros(1,tnum);
fc = 0.05;     % Cutoff frequency (normalised 0.5=nyquist)
n0 = 30;        % Initial delay (samples)
sigma=sqrt(2*log(2))/(2*pi*(fc/dt));
n=0:tnum;
source1=exp(-dt^2*(n-n0).^2/(2*sigma^2)).*(10^-12*10^(80/20));
for n = ceil(tnum/10) : 1 : ceil(tnum/10) + 9 
source1(n) = source1((n-1) * 2);       
source2(n) = source2((n-1) * 2);
end
        
% initialize the velocity and pressure matrices (matrices are set up in a
% y by x fashion to properly display the 2D space (y = rows, x = columns))
p = zeros(ycells - 1, xcells - 1);
ux = zeros(ycells - 1, xcells);
uy = zeros(ycells, xcells - 1);

% set up the multiplication constants for the update equations
uCx = dt/(gx*rho);
uCy = dt/(gy*rho);
pCx = c^2*rho*dt/gx;
pCy = c^2*rho*dt/gy;

% set the characteristic impedances of the walls
ZR = rho*c*(1 + sqrt(1 - alphaR))/(1 - sqrt(1 - alphaR));
ZL = rho*c*(1 + sqrt(1 - alphaL))/(1 - sqrt(1 - alphaL));
ZT = rho*c*(1 + sqrt(1 - alphaF))/(1 - sqrt(1 - alphaF));
ZB = rho*c*(1 + sqrt(1 - alphaB))/(1 - sqrt(1 - alphaB));

% calulcate the coefficients used for the boundary conditions
Rx = rho*gx/dt;
Ry = rho*gy/dt;

% plot vectors
linex = linspace(0, lx - gx, xcells-1);
liney = linspace(0, ly - gx, ycells-1);

%Initialize recording vectors
leftear = zeros(1,tnum);
rightear = zeros(1,tnum);
% loop to update the velocities and pressures over the time steps, n
n = 1;
% while or((max(max(abs(p(:,:)))) > (p0 * 10^(40/10))),(n < 48000))

while n <= (T/dt)

    n = n + 1;
    tic();
    [p, ux, uy] = FDTD2Dfun(p, pCx, pCy, ux, uy, uCx, uCy, Rx, Ry, ZL, ZR, ZT, ZB);
    % set the pressure at the source location
    % NOTE: source vectors for unused drivers will be zeros
    p(s1loc(1),s1loc(2)) = p(s1loc(1),s1loc(2)) - source1(n);
    exectime(n) = toc();
    reciever(n) = p(recieverleftloc(1),recieverleftloc(2));
    
    %PLOTTING SECTION
        surf(linex, liney, p);
        colormap('winter');
        shading interp;
        title(sprintf('Time = %.6f s',n*dt),...
            'Color',[0 0 0],'FontSize', 14);
        xlabel('Width (meters)', 'Color', [0 0 0]);
        ylabel('Length (meters)', 'Color', [0 0 0]);
        caxis([-1e-10 1e-10]);
        drawnow;     
end

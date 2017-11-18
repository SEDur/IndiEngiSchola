
%% Initz Matlab
% clear all;
% %close all;
addpath(genpath('../mls/mls/'));
figure(1)
set(1, 'windowstyle','docked','color', 'w');
figure(2)
set(2, 'windowstyle','docked','color', 'w');

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
fmax = 10000 * hertz;

%grid size
gx = c * (1/fmax) / cstab;
gy = c * (1/fmax) / cstab;

%Dims
%Dim Size (m)
lx = 30*meters;
ly = 30*meters;
xcells = ceil(lx/gx);
ycells = ceil(lx/gy);

%Boundary Absorption Coefs (0 to 1)
alphaL = 1.0;
alphaR = 1.0;
alphaF = 1.0;
alphaB = 1.0;

%number of sources
snum = 2;
%source locations
s1loc = [ceil((11/gx)) ceil((40/gx)/2)];
s2loc = [(s1loc(1) - ceil(1.05/gx)) (s1loc(2) + ceil(1.05/gx))];%source frequency

%recieves position
recieverleftloc = [ceil(ycells/2) ceil(xcells/2)];

%Time of sim
dt = 1/ (c*sqrt((1/(gy^2))+(1/(gx^2))));
T = 1;

% generate the source(s) & determine number of time steps needed
tnum = ceil(T/dt);
source1 = zeros(1,tnum);
source2 = zeros(1,tnum);

stimulus = create_stimulus('tone', 1e-4, 60, 1, 0, [], 1/dt);
gradDelay = ceil(4.17e-3 / dt);
source1 = stimulus;
source2 = [zeros(gradDelay,1); -stimulus];

% initialize the velocity and pressure matrices (matrices are set up in a
% y by x fashion to properly display the 2D space (y = rows, x = columns))
p = zeros(ycells - 1, xcells - 1);
ux = zeros(ycells - 1, xcells);
uy = zeros(ycells, xcells - 1);

pObs = zeros(ycells - 1, xcells - 1);
uxObs = zeros(ycells - 1, xcells);
uyObs = zeros(ycells, xcells - 1);

% set up the multiplication constants for the update equations
uCx = dt/(gx*rho);
uCy = dt/(gy*rho);
pCx = c^2*rho*dt/gx;
pCy = c^2*rho*dt/gy;

% characteristic impedances of the walls
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

n = 1;
% tic
while n <= (T/dt)
    n = n + 1;
    [pObs, uxObs, uyObs] = FDTD2DfunO3(pObs, pCx, pCy, uxObs, uyObs, uCx, uCy, Rx, Ry, ZL, ZR, ZT, ZB, gx);
    [p, ux, uy] = FDTD2Dfun(p, pCx, pCy, ux, uy, uCx, uCy, Rx, Ry, ZL, ZR, ZT, ZB);
    
    p(s1loc(1),s1loc(2)) = p(s1loc(1),s1loc(2)) - source1(n);
    p(s2loc(1),s2loc(2)) = p(s2loc(1),s2loc(2)) - source2(n);
    
    pObs(s1loc(1),s1loc(2)) = pObs(s1loc(1),s1loc(2)) - source1(n);
    pObs(s2loc(1),s2loc(2)) = pObs(s2loc(1),s2loc(2)) - source2(n);
    
    reciever(n) = p(recieverleftloc(1),recieverleftloc(2));
    
    %     rmsP = (rmsOneMinusAlpha .* abs(p) + rmsAlpha .* stateP);
    %     rmsP = 1.0 .* (p) + -0.5 .* stateP;
    %     stateP = rmsP;
    %     rmsP = sqrt(rmsP;
    %     rmsP = sqrt(rmsP);
    
    %PLOTTING SECTION
    figure(1);
    surf(linex, liney, p);
    colormap('winter');
    shading interp;
    title(sprintf('Time = %.6f s',n*dt),...
        'Color',[0 0 0],'FontSize', 14);
    xlabel('Width (meters)', 'Color', [0 0 0]);
    ylabel('Length (meters)', 'Color', [0 0 0]);
    colorbar;
    caxis([-1e-6 1e-6]);
    zlim([-1e-6 1e-6]);
    view(2);
    
    figure(2);
    surf(linex, liney, pObs);
    colormap('winter');
    shading interp;
    title(sprintf('Time = %.6f s',n*dt),...
        'Color',[0 0 0],'FontSize', 14);
    xlabel('Width (meters)', 'Color', [0 0 0]);
    ylabel('Length (meters)', 'Color', [0 0 0]);
    colorbar;
    caxis([-1e-6 1e-6]);
    zlim([-1e-6 1e-6]);
    view(2);
    drawnow;
end

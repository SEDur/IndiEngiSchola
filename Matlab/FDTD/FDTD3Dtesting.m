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
p0 = 10^-12;
rho = 1.21;

%%
%%Hard Code Variables
%Maximum calculation frequency
fmax = 44100 * hertz;

%grid size
gx = c * (1/fmax) / cstab;
gy = c * (1/fmax) / cstab;
gz = c * (1/fmax) / cstab;

%Dims
%Dim Size (m)
lx = 5*meters;
ly = 4*meters;
lz = 3*meters;

xcells = ceil(lx/gx);
ycells = ceil(ly/gy);
zcells = ceil(lz/gz);

%number of sources
snum = 2;
%source locations
sourcelocations = [ceil((ly/gy)/2) ceil((lx/gx)/2) ceil((lz/gz)/2);...
                    ceil((ly/gy)/2)+1 ceil((lx/gx)/2)+1 ceil((lz/gz)/2)+1];
s1Freq = 400;
s2Freq = 400;
%source phase
s1Phase = 0;
s2Phase = 0;
%Source amplitude s
A = 1;

%recieves position
% recieverleftloc = [floor((ycells/2) - (0.1/gy)) ceil((xcells/2)-2) ceil(zcells/2)];
recieverleftloc = [ceil(ycells/2) ceil(xcells/2) ceil(zcells/2)];

recieverrightloc = [ceil((ycells/2) + (0.1/gy)) ceil((xcells/2)+2) ceil(zcells/2)];

%Time of sim
% dt = 1/ (c*sqrt(3/(gx)^2));
dt = 1/ (c*sqrt((1/(gy^2))+(1/(gx^2))+(1/(gz^2))));
% dt = 3.35563e-4;
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
for n = 37 : length(source1)
    if(source1(n) < 0)
       source1(n) = 0; 
    end
end

% source1 = (sin(2*pi*500*[0:dt:T-dt])).*(p0*10^(100/10));
% source2 = (sin(2*pi*500*[0:dt:T-dt])).*(p0*10^(100/10));
        
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
alphaL = 1.0;
alphaR = 1.0;
alphaF = 1.0;
alphaB = 0.1; % This diverges
alphaT = 1.0;
alphaG = 1.0;

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
leftear = zeros(1,tnum);
rightear = zeros(1,tnum);
meanpstore = zeros(1,tnum);

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
while ((n*dt <= T)||(meanpstore(n) > (max(meanpstore)-60)))
    n = n + 1;
    if mod(n,100)
    n*dt
    end    
    [p, ux, uy, uz] = FDTD3Dfun(p, pCx, pCy, pCz, ux, uy, uz, uCx,...
        uCy, uCz, Rx, Ry, Rz, ZL, ZR, ZF, ZB, ZT, ZG);
    p = FDTD3Dsources(p,sourcelocations ,source1(n) , 'soft');
    reciever(n) = p(recieverleftloc(1),recieverleftloc(2),recieverleftloc(3));
%     rightear(n) = p(recieverrightloc(1),recieverrightloc(2),recieverrightloc(3)/p0);
    %PLOTTING SECTION
    signal(n) = real(10*log10(source1(n)/p0));
    meanpstore(n) = 10*log10(mean(mean(mean(abs(real(p)))))/p0);
%     FDTD3Dplotdomain(p, xcells, ycells, zcells, n, dt, p0);
        
end
figure(2);
        subplot(2,1,1);
        plot(dt:dt:n*dt, leftear(1:n));
        hold on;
        plot(dt:dt:n*dt, rightear(1:n));
        plot(dt:dt:n*dt, meanpstore(1:n));
        hold off;
        legend('left','right', 'mean over grid')
        title((sprintf('Current P recieved by listener = %.3f dB & The total sim time was %.6f',(rightear(n)),n*dt)),...
            'Color',[0 0 0],'FontSize', 14);
        ylim([0 max(signal)]);
        subplot(2,1,2);
        plot(dt:dt:n*dt, signal(1:n));
        title('whats sent out by the source');
        ylim([-100 max(signal)]);

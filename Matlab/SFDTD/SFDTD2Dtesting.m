% twoD.m
% S Durbridge 
% 2016

%% Initz Matlab
clear all;
% close all;
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
cstab = 2/(pi*sqrt(2));
%%
%%Hard Code Variables
%Maximum calculation frequency
fmax = 10000 * hertz;
dt = 1/fmax;
%grid size
gx = c * dt / cstab;
gy = c * dt / cstab;

%Dims
%Dim Size (m)
lx = 5*meters;
ly = 4*meters;

pidxRow = [];
pidxCol = [];
uxidx = [];
uyidx = [];

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
% s1loc = [ceil(xcells/3) ceil(ycells/3)];
% s2loc = [ceil(xcells/1.5) ceil(ycells/1.5)];%source frequency
% s1loc = [ceil((lx/gx)/2)-1 ceil((lx/gy)/2)-1];
% s2loc = [ceil((ly/gx)/4) ceil((ly/gy)/2)];%source frequency
sourcelocations = [ceil((1/gy)) ceil(1/gx)];
s1Freq = 400;
s2Freq = 400;
%source phase
s1Phase = 0;
s2Phase = 0;
%Source amplitude 
A = 1;

%recieves position
% recieverleftloc = [ceil(lx/gx/2.42) ceil(ly/gx/8)];
% recieverrightloc = [ceil(lx/gx/2.27) ceil(ly/gx/8)];
recieverleftloc = [ceil(ycells/2) ceil(xcells/2)];

recieverrightloc = [ceil((ycells/2) + (0.1/gy)) ceil((xcells/2)+2)];


%Time of sim
% dt = 1/ (c*sqrt(3/(gx)^2));
% dt = 1/ (c*sqrt((1/(gx^2))+(1/(gy^2))));
% dt = 3.35563e-4;
T = 0.3 ;

% generate the source(s) & determine number of time steps needed

tnum = ceil(T/dt);
source1 = zeros(1,tnum);
source2 = zeros(1,tnum);
% %           t0 = ceil(T/dt) + 1;
%             t0 = 10;
%             t1 = 0 : dt : T;
%             phi = s1Phase*pi;
%             y = A*sin(2*pi*s1Freq*t1 + phi);
%             gain = linspace(0, 1, ceil(length(y)/10));
%             temp = ones(1, length(y));
%             temp(1 : ceil(length(y)/10)) = gain;
%             y = y.*temp;
%             source1(1, t0 : t0 + ceil(T/dt) - 1) = y;
% source1(1,11:1811) = (sin(0:(pi/1800)*2:(2*pi)))*(p0*10^(100/10));
% source2(1,11:1811) = (sin(0:(pi/1800)*2:(2*pi)))*(p0*10^(100/10));
% source1(1,11:20) = ones(1,10).*(10^-12*10^(90/10));
% for n = ceil(tnum/10) : 1 : ceil(tnum/10) + 9 
% source1(n) = source1((n-1) * 2);       
% source2(n) = source2((n-1) * 2);
% end
% 
% source1 = decimate(source1, 2, 'fir')';
% source1 = interp(source1, 2);
% %             t0 = ceil(T/dt) + 1;
%             t0 = 10;
%             t1 = 0 : dt : T;
%             phi = s2Phase*pi;
%             y = A*sin(2*pi*s2Freq*t1 + phi);
%             gain = linspace(0, 1, ceil(length(y)/10));
%             temp = ones(1, length(y));
%             temp(1 : ceil(length(y)/10)) = gain;
%             y = y.*temp;
%             source2(1, t0 : t0 + ceil(T/dt) - 1) = y;

fc = 0.05;     % Cutoff frequency (normalised 0.5=nyquist)
n0 = 30;        % Initial delay (samples)
sigma=sqrt(2*log(2))/(2*pi*(fc/dt));
n=0:tnum;
source1=exp(-dt^2*(n-n0).^2/(2*sigma^2)).*((2*10^-5)*10^(100/20));
source1= source1 ./ max(source1);

% initialize the velocity and pressure matrices (matrices are set up in a
% y by x fashion to properly display the 2D space (y = rows, x = columns))
p = zeros(ycells - 1, xcells - 1);
ux = zeros(ycells - 1, xcells);
uy = zeros(ycells, xcells - 1);

idx3 = uy;

% set up the multiplication constants for the update equations
uCx = dt/(gx*rho);
uCy = dt/(gy*rho);
pCx = c^2*rho*dt/gx;
pCy = c^2*rho*dt/gy;

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
Rx = rho*gx/dt;
Ry = rho*gy/dt;

% plot vectors
linex = linspace(0, lx - gx, xcells-1);
liney = linspace(0, ly - gy, ycells-1);

%Initialize recording vectors
leftear = zeros(1,tnum);
rightear = zeros(1,tnum);

[xvec, yvec] = meshgrid(0 : gx : gx * (xcells-2),...
    0 : gy : gy * (ycells-2));

% loop to update the velocities and pressures over the time steps, n
n = 1;
% while or((max(max(abs(p(:,:)))) > (p0 * 10^(40/10))),(n < 48000))
for n = 1:T/dt    
% n = n + 1;
    tic;
    [idx] = SPARSEfun2DB(p, 0, p0);
    [p, ux, uy] = SFDTD2Dfun(p, pCx, pCy, ux, uy, uCx, uCy, Rx, Ry, ZL,...
        ZR, ZT, ZB, idx);
    extime(n) = toc;
    % set the pressure at the source location
    % NOTE: source vectors for unused drivers will be zeros
    p(sourcelocations(1),sourcelocations(2)) = p(sourcelocations(1),sourcelocations(2)) - source1(n);
%     p(s2loc(1),s2loc(2)) = p(s2loc(1),s2loc(2)) + -source2(n);
%     power(n) = 20*log10(abs(max(p)));
%     leftear(n) = abs(p(recieverleftloc(1),recieverleftloc(2)));
%     rightear(n) = abs(p(recieverrightloc(1),recieverrightloc(2)));
      reciever(n) = p(recieverleftloc(1),recieverleftloc(2));

%PLOTTING SECTION
%         figure(1);
        surf(linex, liney, abs(p));
%         shading interp;
        title(sprintf('Time = %.6f s, Executes at %.6f s',n*dt,extime(n)),...
            'Color',[0 0 0],'FontSize', 14);
        xlabel('Width (meters)', 'Color', [0 0 0]);
        ylabel('Length (meters)', 'Color', [0 0 0]);
%         view(2);
% %         view([25.6 61.2]);
        drawnow;
        
end
% leftear = real(10*log10(leftear/p0));
% rightear = real(10*log10(rightear/p0));
% signal = real(10*log10(source1/p0));
% 
% ax = gca;
% ax = plot(leftear);
% hold on;
% plot(rightear);
% ax.XTickLabel = [0 : (dt)* 10 : length(leftear)* dt];
% hold off;

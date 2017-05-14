% FDTD1D Testing
% S Durbridge 
% 2016

%% Initz Matlab
% clear all;
% close all;
figure(1);
set(1,'color', 'w','windowstyle','docked');

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
p0 = 10^-12;
cstab = sqrt(1/3);
%%
%%Hard Code Variables
%Maximum calculation frequency
fmax = 1000 * hertz;
%grid size
gx = c * (1/fmax) / cstab;
%Dims
%Dim Size (m)
lx = 10*meters;

xcells = ceil(lx/gx);

%Boundary Absorption Coefs (0 to 1)
alphaL = 1.0;
alphaR = 1.0;

%number of sources
snum = 2;
%source locations
s1loc = 2;
s1Freq = 400;
%source phase
s1Phase = 0;
%Source amplitude 
A = 1;

%recieves position
recieverloc = xcells - 1;

% dt = 1/fs;
% %dfine grid spacing
% dx = c * sqrt(2) * dt;
% dx/c * sqrt(2) = dt;


%Time of sim
dt = gx / (c * sqrt(2));
% dt = 1/ (c*sqrt((1/(gx^2))+(1/(gy^2))));
% dt = 3.35563e-4;

T = 1*seconds ;

% generate the source(s) & determine number of time steps needed

tnum = ceil(T/dt);
source1 = zeros(1,tnum);
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
% source1(1,1:1801) = (sin(0:(pi/1800)*2:(2*pi)))*(p0*10^(80/20));
source1(10:1010) = (p0*10^(80/20)) * sin(2*(pi/1010)*(1:1001));
for n = ceil(tnum/10) : 1 : ceil(tnum/10) + 9 
source1(n) = source1((n-1) * 2);
end
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
        
% initialize the velocity and pressure matrices (matrices are set up in a
% y by x fashion to properly display the 2D space (y = rows, x = columns))
p = zeros(1, xcells - 1);
ux = zeros(1, xcells);

% set up the multiplication constants for the update equations
uCx = dt/(gx*rho);
pCx = c^2*rho*dt/gx;
% set the wall reflection coefficients
% if alphaX = 0, then slightly adjust to avoid infinite characteristic
% impedance.
if alphaR == 0
   alphaR = 1e-016; 
end
if alphaL == 0
   alphaL = 1e-016; 
end

% set the characteristic impedances of the walls
ZR = rho*c*(1 + sqrt(1 - alphaR))/(1 - sqrt(1 - alphaR));
ZL = rho*c*(1 + sqrt(1 - alphaL))/(1 - sqrt(1 - alphaL));

% calulcate the coefficients used for the boundary conditions
Rx = rho*gx/dt;

% plot vectors
linex = linspace(0, lx - gx, xcells-1);

%Initialize recording vectors
leftear = zeros(1,tnum);
% loop to update the velocities and pressures over the time steps, n
n = 1;
while or((max(max(abs(p(:,:)))) > (p0 * 10^(40/10))),(n < 48000))
    n = n + 1;
    if mod(n,100)
    (100/tnum)*n;
    10*log10(real(max(abs(p)))/p0)
    end
    
    [p, ux] = FDTD1Dfun(p, pCx, ux, uCx, Rx, ZL, ZR);
    % set the pressure at the source location
    % NOTE: source vectors for unused drivers will be zeros
    p(s1loc) = p(s1loc) - source1(n);
%     p(s2loc(1),s2loc(2)) = p(s2loc(1),s2loc(2)) + -source2(n);
%     power(n) = 20*log10(abs(max(p)));
    leftear(n) = abs(p(recieverloc));
    %PLOTTING SECTION

        plot(linex, abs(p));
        title(sprintf('Time = %.6f s',n*dt),...
            'Color',[0 0 0],'FontSize', 14);
        xlabel('Width (meters)', 'Color', [0 0 0]);
        ylim([0 1e-8])
        drawnow;

end
leftear = real(10*log10(leftear/p0));
signal = real(10*log10(source1/p0));
figure;
ax = gca;
ax = plot(leftear);
hold on;
plot(rightear);
ax.XTickLabel = [0 : (dt)* 10 : length(leftear)* dt];
hold off;

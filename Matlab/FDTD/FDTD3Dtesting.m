%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FDTD3Dtesting.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% A script that will execute an FDTD simulation in 3D.
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initz Matlab
% clear all;
% close all;
%Units
figure(1);
set(1, 'WindowStyle', 'Docked');
cstab = 2/(pi*sqrt(3));
c = 343;
p0 = 2*10^-5;
rho = 1.21;

%% Variable Initialization
%Hard Code Variables
%Maximum calculation frequency
fmax = 2000 * hertz; %% << fs for nyquist would be twice this

%grid discretization step (m)
dx = (c / (6*fmax));
dy = (c / (6*fmax));
dz = (c / (6*fmax));

%timestep
dt = ((1/c)*dx)/2;

%Dims
%Dim Size (m)
lx = 5*meters;
ly = 4*meters;
lz = 3*meters;

%Total number of spatial steps in each direction
xcells = ceil(lx/dx);
ycells = ceil(ly/dy);
zcells = ceil(lz/dz);

%number of sources
snum = 1; 

%source locations
sourcelocations = [ceil((1/dy)) ceil(1/dx) ceil(1/dz);...
                    ceil((1/dy)) ceil(1/dx) ceil(1/dz)];

%recieves position
recieverleftloc = [ceil(ycells/2) ceil(xcells/2) ceil(zcells/2)];

%Time of sim
T = 1.0;

%Generate Source Terms
%Tone Burst

tone = dsp.SineWave('Amplitude',((2*10^-5)*10^(100/20)),...
    'Frequency', 1000,...
    'SampleRate', 1/dt,...
    'SamplesPerFrame',0.01/dt);
w1 = window(@gausswin,0.01/dt,2.5); 
toneBurst = tone() .* w1;
source1 = zeros(T/dt,1);
source1(10:129) = toneBurst;
source1(410:529) = toneBurst;
source1(810:929) = toneBurst;

%MLS
% source1 = GenerateMLSSequence(3,11,0).*((2*10^-5)*10^(100/20));
% T = length(source1)*dt;
% w1 = window(@gausswin,length(source1),2.5); 
% source1 = source1 .* w1;


% pressure and velocity matricies defining the limits of the domain
%(y = rows, x = columns))
p = zeros(ycells - 1, xcells - 1, zcells - 1);
ux = zeros(ycells - 1, xcells, zcells - 1);
uy = zeros(ycells, xcells - 1, zcells - 1);
uz = zeros(ycells - 1, xcells - 1, zcells);

% set up the multiplication constants for the update equations
uCx = dt/(dx*rho);
uCy = dt/(dy*rho);
uCz = dt/(dz*rho);
pCx = c^2*rho*dt/dx;
pCy = c^2*rho*dt/dy;
pCz = c^2*rho*dt/dz;

%Boundary Absorption Coefficients
alphaxn = 0.45;
alphaxp = 0.45;
alphayn = 0.45;
alphayp = 0.45;
alphazp = 0.45;
alphazn = 0.45;

if alphaxp == 0
   alphaxp = 1e-016; 
end
if alphaxn == 0
   alphaxn = 1e-016; 
end
if alphayn == 0
   alphayn = 1e-016; 
end
if alphayp == 0
   alphayp = 1e-016; 
end
if alphazp == 0
   alphazp = 1e-016; 
end
if alphazn == 0
   alphazn = 1e-016; 
end

% set the characteristic impedances of the walls
impedancexn = rho*c*(1 + sqrt(1 - alphaxp))/(1 - sqrt(1 - alphaxp));
impedancexp = rho*c*(1 + sqrt(1 - alphaxn))/(1 - sqrt(1 - alphaxn));
impedanceyn = rho*c*(1 + sqrt(1 - alphayn))/(1 - sqrt(1 - alphayn));
impedanceyp = rho*c*(1 + sqrt(1 - alphayp))/(1 - sqrt(1 - alphayp));
impedancezp = rho*c*(1 + sqrt(1 - alphazp))/(1 - sqrt(1 - alphazp));
impedancezn = rho*c*(1 + sqrt(1 - alphazn))/(1 - sqrt(1 - alphazn));

% calulcate the coefficients used for the boundary conditions
Rx = rho*dx/dt;
Ry = rho*dy/dt;
Rz = rho*dz/dt;

%plot vectors
linex = linspace(0, lx - dx, xcells-1);
liney = linspace(0, ly - dy, ycells-1);
linez = linspace(0, lz - dz, zcells-1);

%plotgrid
[xvec, yvec, zvec] = meshgrid(0 : dx : dx * (xcells-2),...
    0 : dy : dy * (ycells-2), 0 : dz : dz * (zcells-2));
%%
% loop to update the velocities and pressures over the time steps, n
n = 1;
while n*dt < T
    n = n + 1; %Incrememnts n
    
    T - (n * dt) %<< Command Window display of time left to be simulated
    
    tic(); % Begin timer
    
    [p, ux, uy, uz] = FDTD3Dfun(p, pCx, pCy, pCz, ux, uy, uz, uCx,...
        uCy, uCz, Rx, Ry, Rz, impedancexp, impedancexn, impedanceyn,...
        impedanceyp, impedancezp, impedancezn); % Execute Update Equations
    
    p = FDTD3Dsources(p,sourcelocations ,source1(n) , 'soft');
    % Update pressures with the source terms
    
    reciever(n,1) = p(ceil(xcells/4), ceil(ycells/4),ceil(zcells/2));
    reciever(n,2) = p(ceil(xcells/2)+ceil(xcells/4),...
        ceil(ycells/4),ceil(zcells/2));    
    reciever(n,3) = p(ceil(xcells/2), ceil(ycells/2)+ceil(xcells/4),...
        ceil(zcells/2));
    reciever(n,4) = p(ceil(xcells/2+ceil(xcells/4)),...
        ceil(ycells/2)+ceil(xcells/4),ceil(zcells/2));
    reciever(n,5) = p(ceil(xcells/2), ceil(ycells/2),ceil(zcells/2));
    
    srcnorm(n) = p(sourcelocations(1,1),sourcelocations(1,2),...
        sourcelocations(1,3)); % Record source position for normalisation
    
    exectime(n) = toc(); %End timer and record execution time
%     FDTD3Dplotdomain(p, xcells, ycells, zcells, n, dt, p0); 
end

%% Postprocessing

for i = 1 : size(reciever,2)
    lag(i) = getlag(reciever(:,i)',srcnorm);%
    recieverCirc(:,i) = circshift(reciever(:,i)',lag(i));
end

norec = recieverCirc ./ max(abs(recieverCirc));
% recanal = AnalyseMLSSequence(reciever',0,3,11,0,0);
[lpsd, lf] = pwelch(norec,hann(5000),[],5000,1/dt);
srcnrm = srcnorm ./ max(abs(srcnorm));
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
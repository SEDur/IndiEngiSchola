%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PTSD1Dtesting.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% A script that will execute an acoustic PSTD simulation in 1D
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initz
clf;
clc;
% clear all;
% close all;
figure(1);
set(1,'windowstyle','docked');
%% Make Variables

%define FS
fs = 4000.0;
%define density
rho = 1.21;
%define speed of sound
c = 343.0;
%define total time
T = 20.0;
%define grid width
gridWidth = 20.0;
%define timestep
dt = 1/fs;
%dfine grid spacing
dx = c * sqrt(2) * dt;
%calculate pconst
pconst = rho * c^2 * (dt/dx);
%calculate uconst
% uconst = (1/rho)*(dt/dx);
uconst = dt/(dx*rho);

%define pml depth 
PMLdepth = 60;
%calc grid size
N = ceil(abs(gridWidth/dx)+2*PMLdepth);
%calculate differentiation matrix
tempdiffmatrix = zeros(1,N);
temp = zeros(1, N);

%Boundary Absorption Coefs (0 to 1)
alphaL = 0.8;
alphaR = 0.8;

%Calc source
src = zeros(1,ceil(T/(dt/2))+1);
src(10:1010) = (10^-12*10^(50/20)) * sin(2*(pi/1010)*(1:1001));
% tnum = ceil(T/dt);
% fc = 0.1;     % Cutoff frequency (normalised 0.5=nyquist)
% n0 = 120;        % Initial delay (samples)
% sigma=sqrt(2*log(2))/(2*pi*(fc/dt));
% n=0:tnum;
% src=exp(-dt^2*(n-n0).^2/(2*sigma^2)).*((10^-12)*(10^(80/20)));
% for n = 1 : length(src)
%     if(src(n) < 0)
%        src(n) = 0; 
%     end
% end
% src = decimate(src, 2, 'fir');
% src = interp(src, 2);
% srcloc = PMLdepth+1;
% srcloc = ceil(N/2);
srcloc = ceil(N/2);

reciever = zeros(1,ceil(T/dt));
%calculate geometry matricies
phat = zeros(1,N);
uhat = zeros(1,N);
pdiffhat = zeros(1,N);
udiffhat = zeros(1,N);
pd = zeros(1,N);
ud = zeros(1,N);

% pd(PMLdepth:N-PMLdepth) = (10^-12*10^(50/20)) * sin(2*(pi/208)*(1:208));

linex = linspace(0, gridWidth - dx, N);

    for i2 = 1 : N-1
        if i2 <  ceil((N-2)/2)
            tempdiffmatrix(i2) =  (i2-1);
        end
        if i2 ==  ceil((N-1)/2)
            tempdiffmatrix(i2) = 0 * (1+0j);
        end
        if i2 >  ceil((N-1)/2)
            tempdiffmatrix(i2) = (i2 - 1 - N) ;
        end
    end
diffmatrix = 1i * tempdiffmatrix;
PMLconst = ones(1,N);
PMLconst = PMLconst .* (3.142*N);
PMLdiff = zeros(1,N);
PMLdiff(1:PMLdepth) = 1:PMLdepth;
PMLdiff(N-(PMLdepth-1):end) = N -(PMLdepth-1) : N;
PMLdiff(1:PMLdepth) = (1.0/3.0).*(((PMLdepth-PMLdiff(1:PMLdepth))./PMLdepth).^3);
PMLdiff(N-(PMLdepth-1) : end) = (1.0/3.0)*(((PMLdiff(N-(PMLdepth-1) : end)-(N-(PMLdepth-1)))./PMLdepth).^3);
% PMLdiff = zeros(1, N);
PMLalphau = uconst*(1./(1+PMLdiff));
PMLalphap = pconst*(1./(1+PMLdiff));
PMLdiff = ((1.0-PMLdiff)./(1.0+PMLdiff));

S = c*(dt/dx);
Rleft = 1 - alphaL;
Rright = 1 - alphaR;
xiL = (1 + Rleft) / (1 + Rleft - 2 * S * Rleft);
xiR = (1 + Rright)/ (1 + Rright - 2 * S * Rright);
z = 0;


%% solve for some time
% tic();
for i = 1 : T/dt
   [pd, ud] = PSTD1Dfun(pd, ud, diffmatrix,...
     PMLdiff, PMLalphau, PMLalphap, PMLconst);
    [pd, ud] = PTSD1Dboundary(pd, ud, PMLdepth, xiL, xiR);
    pd = PTSD1Dsrc(pd, src(i), srcloc);
    %PLOTTING SECTION
    plot(linex, real(pd));
    ylim([-10e-9 10e-9]);
    title(sprintf('Time = %.6f s',dt*i));
    drawnow;
    reciever(i) = pd(floor(N-PMLdepth)-1);
end
% toc();

%% Display the results


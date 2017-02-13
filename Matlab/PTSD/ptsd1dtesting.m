%% PTSD1D testing script

%% Initz
clc;
clear all;
close all;

%% Make Variables

%define FS
fs = 20000.0;
%define density
rho = 1.21;
%define speed of sound
c = 343.0;
%define total time
T = 10.0;
%define grid width
gridWidth = 100.0;
%define timestep
dt = 1/(2*fs);
%dfine grid spacing
dx = 2 * dt * c;
%calculate pconst
pconst = rho * c^2 * (dt/dx) * dt * c;
%calculate uconst
uconst = (1/rho)*(dt/dx)*dt*c;
%define pml depth 
PMLdepth = 30;
%calc time steps
timestep = abs(T/dt);
%calc grid size
N = ceil(abs(gridWidth/dx)+2*PMLdepth);
%calculate differentiation matrix
tempdiffmatrix = zeros(1,N);
temp = zeros(1, N);
%Calc source
src = zeros(1,ceil(T/dt)+1);
src(10:1010) = 10^-12*50*10^(50/20) * sin(2*(pi/1010)*(1:1001));
srcloc = ceil(N/2);
% alpha = 0;
%calculate geometry matricies
phat = zeros(1,N);
uhat = zeros(1,N);
pdiffhat = zeros(1,N);
udiffhat = zeros(1,N);
pd = zeros(1,N);
ud = zeros(1,N);
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
PMLdiff((N-PMLdepth+1):end) = 1 : PMLdepth;
PMLdiff(1:PMLdepth) = (1/3.0).*(((PMLdepth-PMLdiff(1:PMLdepth))./PMLdepth).^3);
PMLdiff((N-PMLdepth+1):end) = (1/3.0).*(PMLdiff((N-PMLdepth+1):end)./PMLdepth).^3;
PMLalphau = uconst*(1./(1+PMLdiff));
PMLalphap = pconst*(1./(1+PMLdiff));
PMLdiff = ((1-PMLdiff)./(1+PMLdiff));
% PMLdiffu = ((1-PMLdiff)/(1+PMLdiff)) - uconst;
% PMLdiffp = ((1-PMLdiff)/(1+PMLdiff)) - pconst;


%% solve for some time
tic();
for i = 1 : T/dt
   [pd, ud] = PSTD1Dfun(pd, ud, diffmatrix,...
     PMLdiff, PMLalphau, PMLalphap, PMLconst);
    pd = PTSD1Dsrc(pd, src(i), srcloc);
    plot(real(pd));
    title(sprintf('Time = %.6f s',dt*i));
    drawnow;
end
toc();

%% Display the results


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PTSD2Dtestinggpu.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% A script that will execute the PSTD on an NVIDIA CUDA enabled GPU. This
% was not included in the study, but shows just how easy it is to do basic
% GPU computing in matlab.
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initz
clc;
clear all;
% close all;
addpath(genpath('../mls/mls/'));
%% Make Variables
g = gpuDevice(1);

%alpha 
a = 1.0;

%define FS
fs = 10000.0;

%define density
rho = 1.21;

%define speed of sound
c = 343.0;

%define total time
T = 1.0;

%define grid width
gridWidth = 30.0;

%define timestep
dt = (1/fs);

%dfine grid spacing
dx = c * sqrt(2) * dt;

%calculate pconst
pconst = rho * c^2 * (dt/dx);

%calculate uconst
uconst = dt/(dx*rho);

%define pml depth 
PMLdepth = 30;

%calc time steps
timestep = abs(T/dt);

%calc grid size
N = ceil(abs(gridWidth/dx)+2*PMLdepth);

%calculate differentiation matrix
tempdiffmatrix = zeros(1,N);

%Calc source
sStart = 44100 * 40;
src = zeros(1,ceil(T/dt)+1);
src = ((10^-12)*10^(50/10)) .* sin(1000*2*pi*(0:dt:1.0));
win = kaiser(length(src),2.5);
src = src .* win';
clear('win');
src = gpuArray(src);
srcloc = ceil(N/3);

spin = -180 :0.005 : 180;
        
pd = zeros(N,N);
pd = gpuArray(pd);
udx = zeros(N,N);
udx = gpuArray(udx);
udy = zeros(N,N);
udy = gpuArray(udy);

% set up the differentiator
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


[mgx mgy] = meshgrid([tempdiffmatrix tempdiffmatrix]);
[mhx mhy] = meshgrid(-[fliplr(tempdiffmatrix) -fliplr(tempdiffmatrix)]);


diffmatrix =  1i.*(mgx(1:N,1:N)) ;
diffmatrix = gpuArray(diffmatrix);

% set up the PML
PMLconst = ones(N,N);
PMLconst = PMLconst .* (pi*N);
PMLdiff = zeros(N,N);
for i = 1 : N
PMLdiff(i,1:PMLdepth) = 1:PMLdepth;
PMLdiff(i,1:PMLdepth) = (a/3.0).*(((PMLdepth-PMLdiff(i,1:PMLdepth))./PMLdepth).^3);
end
PMLdiff(:,N-PMLdepth+1:end) = fliplr(PMLdiff(:,1:PMLdepth));
PMLdiff2 = PMLdiff';
PMLdiff = sqrt(PMLdiff.^2 + PMLdiff2.^2);
PMLdiffmax = max(max(PMLdiff));
PMLdiffsetmax = 0.3011;
PMLdiff(PMLdiff > PMLdiffsetmax) = PMLdiffsetmax;
PMLalphau = uconst*(1./(1+PMLdiff));
PMLalphau = gpuArray(PMLalphau);
PMLalphap = pconst*(1./(1+PMLdiff));
PMLalphap = gpuArray(PMLalphap);
PMLdiff = ((1-PMLdiff)./(1+PMLdiff));
PMLdiff = gpuArray(PMLdiff);
%% solve for some time
for i = 1 : T/dt
   [pd, udx, udy] = PSTD2Dfun(pd, udx, udy, diffmatrix,...
     PMLdiff, PMLalphau, PMLalphap, PMLconst, N);
    pd = PTSD2Dsrc(pd, src(i), srcloc);
    reciever(i) = abs(gather(pd(ceil(N/2), ceil(N/2))));
    %PLOTTING SECTION
    if mod(i, 100) < 1
    localpd = gather(pd);
    mesh(abs(localpd)); 
    shading interp;
    title(sprintf('Time = %.6f s',dt*i));
    drawnow;
    end
end
reset(g);
plot(reciever);


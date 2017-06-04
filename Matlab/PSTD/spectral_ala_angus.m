
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% spectral_ala_angus.m
% Created by S Durbridge as part of work on a masters dissertation
% Copywrite S Durbridge 2017
%
% A script for 1D PSTD that was adapted to matlab from the C++ code in the
% appendecies of the work by Angus and Caunce on PSTD applied to the GPU.
%
% Any copies of this function distributed by the autor are done so
% without any form of warranty, and should not be reproduced without
% permission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clc;
% clear all;
% close all;
figure(1);
set(1, 'WindowStyle', 'Docked')
%define FS
fs = 2000.0;

%define density
rho = 1.21;

%define speed of sound
c = 343.0;

%define total time
T = 2.0;

%define grid width
gridWidth = 34.3;
dt = 1 / fs;
dx = c * sqrt(2) * dt;

%calculate pconst
pconst = rho * c^2 * (dt/dx);

%calculate uconst
uconst = (1/rho)*(dt/dx);

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
src = zeros(1,ceil(T/(dt/2))+1);
src(10:1010) = (10^-12*10^(50/20)) * sin(2*(pi/1010)*(1:1001));
alpha = 0;

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
cntr = 1;
%calculate propagation
for i = 0 : dt/2 : T
    phat = fft(pd);
    temp = phat .* diffmatrix;
    pdiffhat = ifft(temp);
    for i2 = 1 : length(pdiffhat)
        if i2 < PMLdepth
           alpha = (1/3)*(((PMLdepth-i2)/ PMLdepth)^3); 
        elseif i2 > N - PMLdepth
            alpha = (1/3)*(((i2 -(N-(PMLdepth-1)))/PMLdepth)^3);
        else
            alpha = 0;
        end
        ud(i2) = ud(i2) * ((1-alpha)/(1+alpha))-uconst * (1/(1+alpha))*(pdiffhat(i2)/(3.142*N));
    end
    uhat = fft(ud);

        temp = uhat .* diffmatrix;

    udiffhat = ifft(temp);
    for i2 = 1 : length(udiffhat)
        if i2 < PMLdepth
           alpha = (1/3)*(((PMLdepth-i2)/ PMLdepth)^3); 
        elseif i2 > N - PMLdepth
            alpha = (1/3)*(((i2 -(N-(PMLdepth-1)))/PMLdepth)^3);
        else
            alpha = 0;
        end
        pd(i2) = pd(i2) * ((1-alpha)/(1+alpha))-pconst * (1/(1+alpha))*(udiffhat(i2)/(3.142*N));
    end
    pd(PMLdepth+1) = pd(PMLdepth+1) +  (src(cntr));
    receiver(cntr) = pd(N-PMLdepth-1);
    plot(real(pd));
    title(sprintf('Time = %.6f s',i));
    drawnow();
    cntr = cntr + 1;
end
Nx = 64; % Number of samples collected along first dimension
Ny = 64; % Number of samples collected along second dimension
dx = 1;  % x increment (i.e., Spacing between each column)
dy = .1; % y increment (i.e., Spacing between each row)
x = 0 : dx : (Nx-1)*dx; % x-distance
y = 0 : dy : (Ny-1)*dy; % y-distance
fsig = sin(2*(pi/64)*(1:64));
[dspx dspy] = meshgrid([fsig fsig]);
data_spacedomain = dspx(1:Nx,1:Ny)+ dspy(1:Nx,1:Ny)./(max(max(dspx)).*2); % random 2D matrix
% data_spacedomain = randn(Ny,Nx); % random 2D matrix
Nyq_kx = 1/(2*dx); % Nyquist of data in first dimension
Nyq_ky = 1/(2*dy); % Nyquist of data in second dimension
dkx = 1/(Nx*dx);   % x-Wavenumber increment
dky = 1/(Ny*dy);   % y-Wavenumber increment
kx = -Nyq_kx : dkx : Nyq_kx-dkx; % x-wavenumber
ky = -Nyq_ky : dky : Nyq_ky-dky; % y-wavenumber
% Compute 2D FFT
data_wavenumberdomain = fft2(data_spacedomain);
% Compute grid of wavenumbers
[KX, KY] = meshgrid(ifftshift(kx),ifftshift(ky));
% Compute 2D derivative
data_wavenumberdomain_differentiated = (2i*pi)^2.*KX.*KY.*data_wavenumberdomain; 
% Convert back to space domain
data_spacedomain_differentiated = ifft2(data_wavenumberdomain_differentiated );
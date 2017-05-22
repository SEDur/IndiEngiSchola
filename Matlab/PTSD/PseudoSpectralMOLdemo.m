
 function PseudoSpectralMOLdemo( varargin )
%PseudoSpectralMOLdemo - uses ode45 to integrate wave equations by MOL
%
% PseudoSpectralMOLdemo(m)
%
% solves 1-D wave equation on periodic domain x=[0,1] using ODE45 and
%  sparse matrix representation for RHS -c*D*u
%
% m: number of grid points in domain:  periodicity assumes u(m+1)=u(m)
% if m not passed defaults to 100
% n: number of points in stencil

%set number of grid points
m = 100;
if nargin == 1
    m=varargin{1};
end

% wave speed
c = 343;

% maximum run time (1 period)
n_periods=5.;
tmax = n_periods/c;

% calculate domain
x=linspace(0.,1.,m+1); % m+1 is periodic
x = x(1:end-1); % clip off last point
h = x(2)-x(1) % grid spacing


% initial condition function : gaussian centered at x0
eta = @(x,amp,x0,sigma) amp*exp(-((x-x0)/sigma).^2);
amp = 2.;
sigma = .05;
x0 = .5;
u0 = eta(x,amp,x0,sigma);

% u0 = u0 .* [1:m];

% RHS function
func = @(t,u) -c*psD(u,h);

% set options and run ODE45
options=odeset('RelTol',1.e-4,'AbsTol',1.e-6 ,'InitialStep',.5*h/c);
tic;
[ t, u] = ode45(func,[0 tmax],u0,options);
time = toc;


dt=diff(t);
courant = c*dt/h;
disp(sprintf('Courant #: mean=%g, min=%g, max=%g',mean(courant),min(courant),max(courant)));

% plot out;
prettyplot(t,u,x,eta,amp,x0,sigma,courant);

 
% Calculate error

xmax = 1.;
x0_end = mod(x0+c*t(end),xmax);
u_true = eta(x,amp,x0_end,sigma);
abs_err = h*norm(u(end,:)-u_true,2);
rel_err = abs_err/h/norm(u_true,2);

fprintf('Elapsed Time =%f s, N=%d points, n=%d steps, abs_err=%g, rel_err=%g\n',time,length(x),length(t),abs_err,rel_err);

end
%Utility Functions

% PseudoSpectral 1st derivative operatory
function [Df] = psD(f,h)
% psD - given an array f, with grid-spacing h, calculate the first deriviatve Df using fast-fourier transforms
%  psD(f,h)
%
%  f: regular spaced array

  N=length(f);

  % calculate 1/N times the Nyquist frequency
  kn=2.*pi/N/h;

  % calculate i*k for frequency ordering of fft
  % ( positive wave numbers, then negative wave numbers)

  ik=1i*kn*[0:N/2 -(N/2)+1:-1]';  % calculate ik

  % pseudo-spectral first derivative
  Df=real(ifft(ik.*fft(f)));
end


function prettyplot(t,u,x,eta,amp,x0,sigma,courant)
% prettyplot -  function to make pretty plot of ODE45 output

% plot out courant # (time steps)
figure;
plot(t(2:end),courant);
grid
set(gca,'FontWeight','bold','FontSize',14)
xlabel('time')
ylabel('Courant #')
title('Adaptive Time stepping: ODE45')

% plot out solution and initial condition
xfine = linspace(0,1,301);
ufine = eta(xfine,amp,x0,sigma);
figure;
plot(xfine,ufine,'r');
hold on;
plot(x,u(end,:),'bx-');
grid;
set(gca,'FontWeight','bold','FontSize',14)
legend('u_{true}','u');
xlabel('distance');

% calculate 2-norm of error norm
h=x(2)-x(1);
utrue=eta(x,amp,x0,sigma);
abs_err = h*norm(u(end,:) - utrue,2);
rel_err = abs_err/h/norm(utrue,2);
title(sprintf('PS-MOL demo: ODE45, CD, m=%d , N=%d,\\sigma=%g,||e||_2=%g',length(x),length(t),sigma,rel_err));




end
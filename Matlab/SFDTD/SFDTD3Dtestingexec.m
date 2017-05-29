%3D implementation
% S Durbridge 
% 2016

%%Initz Matlab
clear all;
% close all;
% figure(1);
% figure( 'color', 'w');
% figure(2);
% figure( 'color', 'w');
%%Initz Variables

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
cstab = 2/(pi*sqrt(3));
%% Initz Winodws
figure(1);
set(1, 'WindowStyle', 'Docked')
figure(2);
set(2, 'WindowStyle', 'Docked')
%% Hard Coded Variables

%Maximum calculation frequency
fmax = 500 * hertz;
% dt = 1/fmax;
% %grid size
% gx = c * dt / cstab;
% gy = c * dt / cstab;
% gz = c * dt / cstab;
gx = (c / (6*fmax));
gy = (c / (6*fmax));
gz = (c / (6*fmax));
dt = ((1/c)*gx)/2;
%Dims
% %Dim Size (m)
% lx = 5*meters;
% ly = 4*meters;
% lz = 3*meters;


%Boundary Absorption Coefs (0 to 1)
alphaL = 0.45;
alphaR = 0.45;
alphaF = 0.45;
alphaB = 0.45;
alphaT = 0.45;
alphaG = 0.45;

%number of sources
snum = 2;
%source locations
% s1loc = [ceil(xcells/3) ceil(ycells/3)];
% s2loc = [ceil(xcells/1.5) ceil(ycells/1.5)];%source frequency
% s1loc = [ceil((ly/gy)/2) ceil((lx/gx)/2) ceil((lz/gz)/2)] ;
% s2loc = [ceil((ly/gy)/4) ceil((lx/gx)/2) ceil((lz/gz)/4)];%source frequency
% sourcelocations = [ceil((1/gy)) ceil(1/gx) ceil(1/gz);...
%                     ceil((1/gy)) ceil(1/gx) ceil(1/gz)];

% s1Freq = 400;
% s2Freq = 400;
% %source phase
% s1Phase = 0;
% s2Phase = 0;
% %Source amplitude 
% A = 1;

%recieves position
% recieverleftloc = [ceil(lx/gx/2.42) ceil(ly/gx/8) ceil((lz/gz)/2)];
% recieverrightloc = [ceil(lx/gx/2.27) ceil(ly/gx/8) ceil((lz/gz)/2)];
% recieverleftloc = [floor((ycells/2) - (0.1/gy)) ceil((xcells/2)-2) ceil(zcells/2)];
% recieverrightloc = [ceil((ycells/2) + (0.1/gy)) ceil((xcells/2)+2) ceil(zcells/2)];


%Time of sim
% dt = 1/ (c*sqrt(3/(gx)^2));
% dt = 1/ (c*sqrt((1/(gy^2))+(1/(gx^2))+(1/(gz^2))));
% dt = 3.35563e-4;
T = 1.0 ;

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
%% Source 1 Gauss
% fc = 0.05;     % Cutoff frequency (normalised 0.5=nyquist)
% n0 = 30;        % Initial delay (samples)
% sigma=sqrt(2*log(2))/(2*pi*(fc/dt));
% n=0:tnum;
% source1=exp(-dt^2*(n-n0).^2/(2*sigma^2)).*(10^-12*10^(80/20));
% for n = 37 : length(source1)
%     if(source1(n) < 0)
%        source1(n) = 0; 
%     end
% end
%% Source 1 step
% source1(1,5:7) = ones(1,3).*(10^-12*10^(40/20));
% for n = ceil(tnum/10) : 1 : ceil(tnum/10) + 9 
% source1(n) = source1((n-1) * 2);       
% source2(n) = source2((n-1) * 2);
% end

% source1 = decimate(source1, 2, 'fir')';
% source1 = interp(source1, 2);
% source1=filter([1 -1], [1 -0.995], source1); % DC Blocker
% source1 = source1 .* (p0*10^(100/10));
% for n = ceil(tnum/10) : 1 : ceil(tnum/10) + 9 
% source1(n) = source1((n-1) * 2);       
% source2(n) = source2((n-1) * 2);
% end
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
% figure(4);
% fc = 0.05;     % Cutoff frequency (normalised 0.5=nyquist)
% n0 = 10;        % Initial delay (samples)
% sigma=sqrt(2*log(2))/(2*pi*(fc/dt));
% n=0:tnum;
% source1=exp(-dt^2*(n-n0).^2/(2*sigma^2));
% source1 = -(((2*10^-5)*10^(100/20))*source1);
% source1(source1 > (2*10^-5)) = (2*10^-5);
% plot(0:dt:(length(source1)-1)*dt,source1);
% xlim([0 0.03]);
% drawnow;
% figure(2);
% tone = dsp.SineWave('Amplitude',((2*10^-5)*10^(100/20)),...
%     'Frequency', 1000,...
%     'SampleRate', 1/dt,...
%     'SamplesPerFrame',0.01/dt);
% w1 = window(@gausswin,0.01/dt,2.5); 
% toneBurst = tone() .* w1;
% source1 = zeros(T/dt,1);
% source1(10:129) = toneBurst;
% source1(410:529) = toneBurst;
% source1(810:929) = toneBurst;
source1 = GenerateMLSSequence(2,9,0).*((2*10^-5)*10^(100/20));
T = length(source1)*dt;
w1 = window(@gausswin,length(source1),2.5); 
source1 = source1 .* w1;

%%
widths = [5 10 20 40 60];
for cntr = 1 : 5
    
    lx = widths(cntr);
    ly = lx;
    lz = ly;

xcells = ceil(lx/gx);
ycells = ceil(ly/gy);
zcells = ceil(lz/gz);


sourcelocations = [ceil(1/gy) ceil(1/gx) ceil(1/gz);...
                    ceil(xcells/2) ceil(ycells/2) ceil(zcells/2)];
recieverleftloc = [ceil(ycells/2) ceil(xcells/2) ceil(zcells/2)];

recieverrightloc = [ceil((ycells/2) + (0.1/gy)) ceil((xcells/2)+2) ceil(zcells/2)];

% initialize the velocity and pressure matrices (matrices are set up in a
% y by x fashion to properly display the 2D space (y = rows, x = columns))
p = ones(ycells - 1, xcells - 1, zcells - 1)*(2*10^-5);
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

%Set zsliceloc
[xvec, yvec, zvec] = meshgrid(0 : gx : gx * (xcells-2),...
    0 : gy : gy * (ycells-2), 0 : gz : gz * (zcells-2));

% loop to update the velocities and pressures over the time steps, n
n = 1;
% while or((mean(mean(mean(abs(real(p(:,:,:)))))) > (p0 * 10^(60/10))),(n < (1600)))
% while or((meanpstore(n) > (max(meanpstore)-60)),(n < (1600)))
for n = 1:T/dt    
%     n = n + 1;
%     if mod(n,100)
%     (100/tnum)*n
%     end
tic;
    disp(strcat('size',num2str(lx),'time',num2str(T-(n*dt))));
    [idx] = SPARSEfun3DC(p, 10, p0);
    [p, ux, uy, uz] = SFDTD3DfunC(p, pCx, pCy, pCz, ux, uy, uz, uCx, uCy, uCz, Rx, Ry, Rz, ZL, ZR, ZT, ZB, ZF, ZB, idx);
    reciever(n,1) = p(ceil(xcells/4), ceil(ycells/4),ceil(zcells/2));
    reciever(n,2) = p(ceil(xcells/2)+ceil(xcells/4), ceil(ycells/4),ceil(zcells/2));    
    reciever(n,3) = p(ceil(xcells/2), ceil(ycells/2)+ceil(xcells/4),ceil(zcells/2));
    reciever(n,4) = p(ceil(xcells/2+ceil(xcells/4)), ceil(ycells/2)+ceil(xcells/4),ceil(zcells/2));
    reciever(n,5) = p(ceil(xcells/2), ceil(ycells/2),ceil(zcells/2));
    srcnorm(n) = p(sourcelocations(1,1),sourcelocations(1,2),sourcelocations(1,3));
    % Input source
    p(sourcelocations(1),sourcelocations(2),sourcelocations(3)) = p(sourcelocations(1),sourcelocations(2),sourcelocations(3)) - source1(n);
     exectime(n) = toc();
    %     p(s2loc(1),s2loc(2)) = p(s2loc(1),s2loc(2)) + -source2(n);
%     power(n) = 20*log10(abs(max(p)));
%     leftear(n) = 10*log10(abs(real(p(recieverleftloc(1),recieverleftloc(2),recieverleftloc(3))))/p0);
%     rightear(n) = 10*log10(abs(real(p(recieverrightloc(1),recieverrightloc(2),recieverrightloc(3))))/p0);
%     %PLOTTING SECTION
%     leftear = real(10*log10(leftear/p0));
%     rightear = real(10*log10(rightear/p0));
%     reciever(n) = p(recieverleftloc(1),recieverleftloc(2),recieverleftloc(3));
%     signal(n) = real(10*log10(source1(n)/p0));
%     meanpstore(n) = 10*log10(mean(mean(mean(abs(real(p)))))/p0);
%     SFDTD3Dplotdomain(p, xcells, ycells, zcells, n, dt, p0);
%     patch(isocaps(idx,1),...
%    'FaceColor','interp','EdgeColor','none');
%     p1 = patch(isosurface(idx,1),...
%    'FaceColor','blue','EdgeColor','none');
%     isonormals(idx,p1);
%     view(3); 
%     axis vis3d tight
%     camlight left
%     colormap('jet');
%     lighting gouraud
%     drawnow;
    %         xslice = ceil(xcells/2); 
%         yslice = ycells-1; 
% %         zslice = (ceil(zcells-1)/2);
% %         slice(abs(p),xslice,yslice,zslice)
%         zslice = (ceil(zcells-1)/2);
% %         figure(1);
%         slice(abs(p),[xcells-1 ceil(xcells/2)],ycells-1,zslice) 
% %         shading interp;
%         title(sprintf('Time = %.6f s Max P = %.3f dB',n*dt,10*log10(real(max(max(max(abs(p(:,:,:))))))/p0)),...
%             'Color',[0 0 0],'FontSize', 14);
%         xlabel('Width (meters)', 'Color', [0 0 0]);
%         ylabel('Length (meters)', 'Color', [0 0 0]);
% 
%         drawnow();
end

norec = reciever ./ max(abs(reciever));
recanal = AnalyseMLSSequence(reciever(:,1),0,2,9,0,0);
% norec = Hd(norec);
% [lpsd, lf] = pwelch(norec,hann(5000),[],5000,fs);
[lpsd, lf] = pwelch(norec,hann(1000),[],1000,1/dt);
% clear('Hd');
% Hd = postprocessingDCfilter;
srcnrm = srcnorm ./ max(abs(srcnorm));
% srcnrm = Hd(srcnrm);
% [spsd, sf] = pwelch(srcnrm,hann(5000),[],5000,fs);
[spsd, sf] = pwelch(srcnrm,hann(1000),[],1000,1/dt);
filename = strcat('xwidth',num2str(lx),'.mat');
save(filename,'exectime', 'norec', 'recanal', 'lpsd', 'lf', 'srcnrm', 'spsd', 'sf');
end
%% Some really minor postprocessing

for i = 1 : size(reciever,2)
    lag(i) = getlag(reciever(:,i)',srcnorm);
    recieverCirc(:,i) = circshift(reciever(:,i)',lag(i));
end

% Hd = postprocessingDCfilter;
norec = recieverCirc ./ max(abs(recieverCirc));
% recanal = AnalyseMLSSequence(reciever',0,3,11,0,0);
% norec = Hd(norec);
% [lpsd, lf] = pwelch(norec,hann(5000),[],5000,fs);
[lpsd, lf] = pwelch(norec,hann(5000),[],5000,1/dt);
% clear('Hd');
% Hd = postprocessingDCfilter;
srcnrm = srcnorm ./ max(abs(srcnorm));
% srcnrm = Hd(srcnrm);
% [spsd, sf] = pwelch(srcnrm,hann(5000),[],5000,fs);
[spsd, sf] = pwelch(srcnrm,hann(5000),[],5000,1/dt);
% figure(3);
%         subplot(2,1,1);
%         plot(dt:dt:n*dt, leftear(1:n));
%         hold on;
%         plot(dt:dt:n*dt, rightear(1:n));
%         plot(dt:dt:n*dt, meanpstore(1:n));
%         hold off;
%         legend('left','right', 'mean over grid')
%         title((sprintf('Current P recieved by listener = %.3f dB & The total sim time was %.6f',(rightear(n)),n*dt)),...
%             'Color',[0 0 0],'FontSize', 14);
%         ylim([0 max(signal)]);
%         subplot(2,1,2);
%         plot(dt:dt:n*dt, signal(1:n));
%         title('whats sent out by the source');
%         ylim([-100 max(signal)]);
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

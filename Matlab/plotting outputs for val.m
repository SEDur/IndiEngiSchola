%% Display the results
subplot(4,1,1);
plot(0:dt:((length(reciever)-1)*dt),reciever)
hold on;
plot(0:dt:((length(reciever)-1)*dt),source1(1:length(reciever)))
hold off;
axis('tight')
legend('reciever','source');
title('raw input and output');
subplot(4,1,2);
plot(0:dt:((length(norec)-1)*dt),norec,'--','linewidth',2.0)
hold on;
plot(0:dt:((length(srcnrm)-1)*dt),srcnrm)
hold off;
axis('tight')
legend('reciever','source');
title('normalised input and output');
subplot(4,1,3);
plot(thrlf, db(lpsd),'--','Linewidth',2.0);
hold on;
plot(sf, db(spsd));
hold off;
legend('reciever','source');
grid('on');
title('power spectral density of input and output');
subplot(4,1,4);
plot(0:dt:((length(reciever)-1)*dt),exectime)
axis('tight')
ttlstr = sprintf('computation time per cycle, total time is %d',sum(exectime));
title(ttlstr);

%% Load Files
thr10 = load('sfdtd10thresh.mat');
thr20 = load('sfdtd20thresh.mat');
thr30 = load('sfdtd30thresh.mat');
thr40 = load('sfdtd40thresh.mat');
thr50 = load('sfdtd50thresh.mat');
thr60 = load('sfdtd60thresh.mat');
thr70 = load('sfdtd70thresh.mat');
thr80 = load('sfdtd80thresh.mat');
%% Plot execution time
smoothingcoeff = ones(1, length(thr10.exectime))/length(thr10.exectime);
% filter(smoothingcoeff, 1, tempC);
plot(0:thr10.dt:((length(thr10.reciever)-1)*thr10.dt),filter(smoothingcoeff, 1, thr10.exectime));
hold on;
plot(0:thr10.dt:((length(thr10.reciever)-1)*thr10.dt),filter(smoothingcoeff, 1, thr20.exectime));
plot(0:thr10.dt:((length(thr10.reciever)-1)*thr10.dt),filter(smoothingcoeff, 1, thr30.exectime));
plot(0:thr10.dt:((length(thr10.reciever)-1)*thr10.dt),filter(smoothingcoeff, 1, thr40.exectime));
plot(0:thr10.dt:((length(thr10.reciever)-1)*thr10.dt),filter(smoothingcoeff, 1, thr50.exectime));
plot(0:thr10.dt:((length(thr10.reciever)-1)*thr10.dt),filter(smoothingcoeff, 1, thr60.exectime));
plot(0:thr10.dt:((length(thr10.reciever)-1)*thr10.dt),filter(smoothingcoeff, 1, thr70.exectime));
plot(0:thr10.dt:((length(thr10.reciever)-1)*thr10.dt),filter(smoothingcoeff, 1, thr80.exectime));
hold off;
grid('on');
title('Execution Time for an SFDTD Simulation with Different Thresholds');
legend('10dB','20dB','30dB','40dB','50dB','60dB','70dB','80dB','Location','southeast');
xlabel('Time in Simulation (s)');
ylabel('Execution Time(s)');
set(gcf,'color','w');
%% 
subplot(3,1,1);
plot(thr20.lf, db(thr20.lpsd),'--','Linewidth',2.0);
hold on;
plot(thr20.sf, db(thr20.spsd));
hold off;
legend('reciever','source');
grid('on');
title('Power Spectral Density of Input and Output for 20dB SFDTD Threshold');
xlim([0 1000]);
subplot(3,1,2);
plot(thr40.lf, db(thr40.lpsd),'--','Linewidth',2.0);
hold on;
plot(thr40.sf, db(thr40.spsd));
hold off;
legend('reciever','source');
grid('on');
title('Power Spectral Density of Input and Output for 40dB SFDTD Threshold');
xlim([0 1000]);
subplot(3,1,3);
plot(thr60.lf, db(thr60.lpsd),'--','Linewidth',2.0);
hold on;
plot(thr60.sf, db(thr60.spsd));
hold off;
legend('reciever','source');
grid('on');
title('Power Spectral Density of Input and Output for 60dB SFDTD Threshold');
xlim([0 1000]);


fdx60data.exectime()

%% Load FDTD data
fdx5data = load('xwidth5.mat');
fdx10data = load('xwidth10.mat');
fdx20data = load('xwidth20.mat');
fdx40data = load('xwidth40.mat');
fdx60data = load('xwidth60.mat');

smoothingcoeff2 = ones(1, length(fdx40data.exectime))/length(fdx40data.exectime);

semilogy(0:fdtddt:((length(fdx5data.norec)-1)*fdtddt),filter(smoothingcoeff2, 1, fdx5data.exectime));
hold on;
semilogy(0:fdtddt:((length(fdx5data.norec)-1)*fdtddt),filter(smoothingcoeff2, 1, fdx10data.exectime));
semilogy(0:fdtddt:((length(fdx5data.norec)-1)*fdtddt),filter(smoothingcoeff2, 1, fdx20data.exectime));
semilogy(0:fdtddt:((length(fdx5data.norec)-1)*fdtddt),filter(smoothingcoeff2, 1, fdx40data.exectime));
semilogy(0:fdtddt:((length(fdx5data.norec)-1)*fdtddt),filter(smoothingcoeff2, 1, fdx60data.exectime));

hold off;
grid('on');
title('Execution Time for FDTD Simulations of Different sizes');
legend('5m^3','10m^3','20m^3','40m^3','60m^3','Location','southeast');
xlabel('Time in Simulation (s)');
ylabel('Execution Time(s)');
set(gcf,'color','w');

%%
semilogy(0:fdtddt:((length(fdx5data.norec)-1)*fdtddt),fdx5data.exectime);
hold on;
semilogy(0:fdtddt:((length(fdx5data.norec)-1)*fdtddt), fdx10data.exectime);
semilogy(0:fdtddt:((length(fdx5data.norec)-1)*fdtddt), fdx20data.exectime);
semilogy(0:fdtddt:((length(fdx5data.norec)-1)*fdtddt), fdx40data.exectime);
semilogy(0:fdtddt:((length(fdx5data.norec)-1)*fdtddt), fdx60data.exectime);

hold off;
grid('on');
title('Execution Time for FDTD Simulations of Different sizes');
legend('5m^3','10m^3','20m^3','40m^3','60m^3','Location','southeast');
xlabel('Time in Simulation (s)');
ylabel('Execution Time(s)');
set(gcf,'color','w');


%%

psx5data = load('xwidth5.mat');
psx10data = load('xwidth10.mat');
psx20data = load('xwidth20.mat');
psx40data = load('xwidth40.mat');
psx60data = load('xwidth60.mat');
smoothingcoeff3 = ones(1, length(psx40data.roundtime))/length(psx40data.roundtime);

plot(0:pstddt:((length(psx5data.norec)-1)*pstddt),filter(smoothingcoeff3, 1, psx5data.roundtime));
hold on;
plot(0:pstddt:((length(psx5data.norec)-1)*pstddt),filter(smoothingcoeff3, 1, psx10data.roundtime));
plot(0:pstddt:((length(psx5data.norec)-1)*pstddt),filter(smoothingcoeff3, 1, psx20data.roundtime));
plot(0:pstddt:((length(psx5data.norec)-1)*pstddt),filter(smoothingcoeff3, 1, psx40data.roundtime));
plot(0:pstddt:((length(psx5data.norec)-1)*pstddt),filter(smoothingcoeff3, 1, psx60data.roundtime));

hold off;
grid('on');
title('Execution Time for PSTD Simulations of Different sizes');
legend('5m^3','10m^3','20m^3','40m^3','60m^3','Location','southeast');
xlabel('Time in Simulation (s)');
ylabel('Execution Time(s)');
set(gcf,'color','w');

%%

psx5data = load('xwidth5.mat');
psx10data = load('xwidth10.mat');
psx20data = load('xwidth20.mat');
psx40data = load('xwidth40.mat');
psx60data = load('xwidth60.mat');
smoothingcoeff3 = ones(1, length(psx40data.roundtime))/length(psx40data.roundtime);

plot(0:pstddt:((length(psx5data.norec)-1)*pstddt),psx5data.roundtime);
hold on;
plot(0:pstddt:((length(psx5data.norec)-1)*pstddt), psx10data.roundtime);
plot(0:pstddt:((length(psx5data.norec)-1)*pstddt),psx20data.roundtime);
plot(0:pstddt:((length(psx5data.norec)-1)*pstddt),psx40data.roundtime);
plot(0:pstddt:((length(psx5data.norec)-1)*pstddt),psx60data.roundtime);

hold off;
grid('on');
title('Execution Time for PSTD Simulations of Different sizes');
legend('5m^3','10m^3','20m^3','40m^3','60m^3','Location','northwest');
xlabel('Time in Simulation (s)');
ylabel('Execution Time(s)');
set(gcf,'color','w');
axis('tight')

%%
sfdx5data = load('xwidth5.mat');
sfdx10data = load('xwidth10.mat');
sfdx20data = load('xwidth20.mat');
sfdx40data = load('xwidth40.mat');
sfdx60data = load('xwidth60.mat');
smoothingcoeff3 = ones(1, length(sfdx40data.exectime))/length(sfdx40data.exectime);

plot(0:sfdtddt:((length(sfdx5data.norec)-1)*sfdtddt),sfdx5data.exectime);
hold on;
plot(0:sfdtddt:((length(sfdx5data.norec)-1)*sfdtddt), sfdx10data.exectime);
plot(0:sfdtddt:((length(sfdx5data.norec)-1)*sfdtddt),sfdx20data.exectime);
plot(0:sfdtddt:((length(sfdx5data.norec)-1)*sfdtddt),sfdx40data.exectime);
plot(0:sfdtddt:((length(sfdx5data.norec)-1)*sfdtddt),sfdx60data.exectime);

hold off;
grid('on');
title('Execution Time for SFDTD Simulations of Different sizes');
legend('5m^2','10m^2','20m^2','40m^2','60m^2','Location','northwest');
xlabel('Time in Simulation (s)');
ylabel('Execution Time(s)');
set(gcf,'color','w');
axis('tight')


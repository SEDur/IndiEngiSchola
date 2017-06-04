RT = [1.6 : 0.1 : 10];
% V = [2250: (800000 - 2250)/104 :800000];
v = [30 : ((150-30)/length(RT)) : 150 -((150-30)/length(RT))];
V = v.^3;
scf = [];
for i = 1 : length(RT)
    for i1 = 1 : length(V)
        scf(i,i1) = 2000*sqrt(RT(i)./V(i1));
    end
end
% surf(RT,v,scf)
surf(V,RT,real(scf));
xscale('log');
ylabel('RT60 (s)');
xlabel('V (m^3)');
zlabel('Schroeder Frequency');
% xlim([30 300])
grid('on');
axis('tight');
shading('interp');
view(2);
c = colorbar;
c.Label.String = 'Frequench (Hz)';
title('Schroeder Frequency as a Function of Room Size & RT60');
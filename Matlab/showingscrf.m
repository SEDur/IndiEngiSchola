RT = [0.6 : 0.1 : 11];
V = [2250: (800000 - 2250)/104 :800000];

for i = 1 : 105
for i1 = 1 : 105
scf(i,i1) = 2000*sqrt(RT(i)./V(i1));
end
end
surf(RT,V,scf)
xlabel('RT');
ylabel('V');
zlabel('Schroeder Frequency');
grid('on');
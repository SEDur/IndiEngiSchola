
% Setup Room
c = 343;
% Dims
l = 5;
w = 4;
h = 3;
% Surfaces
xn = l*h;
xp = xn;
yn = w*h;
yp = yn;
zn = l*w;
zp = zn;
S = (xn+yn+zn)*2;
xna = 0.45;
xpa = 0.45;
yna = 0.45;
ypa = 0.45;
zna = 0.45;
zpa = 0.45;
% volume
v = l * w * h;
% AcousticProperties
a = ((xn*xna)+(xp*xpa)+(yn*yna)+(yp*ypa)+(zn*zna)+(zp*zpa))/S;
Sa = S * a;
RT60 = (0.161*v)/(Sa);
nWallHits = 6*log(10)*1/a;
MFP = 4*v/S;
timePerMFP = MFP/c;

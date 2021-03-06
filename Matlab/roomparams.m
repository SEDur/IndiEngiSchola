
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
RT60 = -((0.161*v)/(S*log(1-a)));
nWallHits = 6*log(10)*1/a;
MFP = 4*v/S;
timePerMFP = MFP/c;
fschroder = 2000*sqrt(RT60/v);

% AXIAL MODES
n = 1:9;
fl = c/2 .* n/l;
fw = c/2 .* n/w;
fh = c/2 .* n/h;
% Tangential modes
nl = 1:10;
nw = 1:10;
nh = 1:10;
flw = c/2 * sqrt(bsxfun(@plus, (nl.'/l).^2, (nw/w).^2));
fwh = c/2 * sqrt(bsxfun(@plus, (nw.'/w).^2, (nh/h).^2));
flh = c/2 * sqrt(bsxfun(@plus, (nl.'/l).^2, (nh/h).^2));
% Oblique Modes
flwh = c/2 * sqrt(bsxfun(@plus,(flw(1:2,1:2) * 2/c).^2, permute(nh(1:2)/h,[3 1 2]).^2));

[srt, idx] = sort([fl(:);fw(:);fh(:);flw(:);fwh(:);flh(:);flwh(:)]);

freqs = 0 : 800;
vals = zeros(length(freqs),1);
for i = 1 : length(srt)
    vals(round(srt(i))) = vals(round(srt(i))) + 1;
end
fscr = zeros(length(vals),1);
fscr(ceil(fschroder)) = 2;
% vals(vals>1) = 1;
stem(fscr,'r')
hold('on');
stem(vals);
hold('off');
xlim([0 ceil(fschroder)+20])
title('Theoretical Modes for a 5m by 4m by 3m Room');
xlabel('Frequency (Hz)');
ylabel('Number of Modes');
legend('Schroeder Frequency','Marker','diamond','Modes','Location','northwest');
% vals = vals ./ max(vals);
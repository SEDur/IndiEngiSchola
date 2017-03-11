% prototyping.m

myfun = -100 : 10 : 100;
myfuninterp = interp(myfun, 2);

threshold = 50;

myfuninterpfft = fft(myfuninterp);

myfuninterpfft = myfuninterpfft .* [zeros(1,21) ones(1,21)];

myfuninterpifft = ifft(myfuninterpfft);
plot(myfuninterp)
hold on;
plot(real(myfuninterpifft))
hold off;
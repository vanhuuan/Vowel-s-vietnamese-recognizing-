%[x,fs] = wavread('ederwander_IN_250Hz.wav');

CorrectFactor = 0.986;
threshold = 0.2;

%F0 start test
f  = 250;
fs = 44100;

signal= 0.9*sin(2*pi*f/fs*(0:9999)); 
x=signal';

framed = x(1:2048);

windowed = framed .* hann(length(framed));

FFT = fft(windowed, 2048);

FFT = FFT(1 : size(FFT,1) / 2);

FFT = abs(FFT);

hps1 = downsample(FFT,1);
hps2 = downsample(FFT,2);
hps3 = downsample(FFT,3);
hps4 = downsample(FFT,4);
hps5 = downsample(FFT,5);

y = [];

for i=1:length(hps5)
      Product =   hps1(i)  * hps2(i) * hps3(i) * hps4(i) * hps5(i);
      y(i) = [Product];
end

[m,n]=findpeaks(y, 'SORTSTR', 'descend');

Maximum = n(1);

 %try fix octave error
if (y(n(1)) * 0.5) > (y(n(2))) %& ( ( m(2) / m(1) ) > threshold )

    Maximum  = n(length(n));

end

F0 =  ( (Maximum / 2048) * fs ) * CorrectFactor 

plot(y)
%%
t = 0:0.001:1-0.00;
x = cos(2*pi*100*t)+randn(size(t));
winvec = hamming(length(x)); % hann(length(x));
xdft = fft(x'.*winvec);
plot(abs(xdft))
%%
y = [1,2,3,4,5,6,7,8,9];
x = [];
x = horzcat(x,y(1:3));
x = horzcat(x,y(4:6));
x = horzcat(x,y(7:9))

reshape(x, 3, [])
%%
[y,Fs] = audioread("02FVA.wav");
y= y(1.12*Fs:1.14*Fs);
nfft = 2048;
window = hamming(length(y)); %Hamming window function
FFT = abs(fft(y.*window,nfft));
FFT = FFT(1:length(FFT)/2);
Fa = [0:Fs/nfft:Fs/2-Fs/nfft]; % frequency axis
[peaks, locs] = findpeaks(FFT,'MinPeakDistance',5)
i=1;
F0 = Fa(locs(i+1)) - Fa(locs(i))
figure(1);
plot(FFT(1:length(FFT)/2));
%%
[y,Fs] = audioread("01MDA.wav");
HPS(y(0.55*Fs:0.57*Fs),Fs);
disp(Fs);
function [F0_HPS] = HPS(Frame,Fs)
    window = hamming(length(Frame)); %Hamming window function
    n = 2048 % So diem lay FFT
    FFT = abs(fft(Frame.*window,n));
    FFT1 = downsample(FFT,1);
    FFT2 = downsample(FFT,2);
    FFT3 = downsample(FFT,3);
    FFT4 = downsample(FFT,4);
    FFT5 = downsample(FFT,5);
    hps = [];
    for i=1:length(FFT5)
          Product =   FFT1(i)  * FFT2(i) * FFT3(i) * FFT4(i) * FFT5(i);
          hps(i) = Product;
    end
    [m,n]=findpeaks(hps, 'SORTSTR', 'descend');
    Maximum = n(1);
    F0 =   (Maximum / 2048) * Fs 
    figure(1);
    plot(FFT(1:length(FFT)/2));
end
%%
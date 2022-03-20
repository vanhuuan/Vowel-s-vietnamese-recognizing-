%% read record signal
file = '44.1khz.wav';
[x, Fs] = audioread(file);
%% play audio to speaker
%sound(x,Fs);
%% plot signal
plot(x);
title("bai3-1");
xlabel("time");
ylabel("signal");
%% energy = sum |x[n]*x[n]| 
E = sum(x.^2);
disp("energy = "+ E);
%% power = energy / length
L=length(x);
P = E/L;
disp("Power = "+ P);
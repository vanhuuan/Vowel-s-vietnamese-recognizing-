file = '44.1khz.wav'
[y,Fs] = audioread(file);
s = length(y)/Fs; 
disp("length = " + s+"s");
disp("consist of " +length(y)+" sample");
time=(1:length(y))/Fs;
tiledlayout(2,1)
nexttile
plot(time,y)
title('Time')
nexttile
plot(y)
title('Sample')
sound(y,Fs);
pause(s);
sound(y,Fs*2);
pause(s/2);
sound(y,Fs/2);
pause(s*2);
disp('done');



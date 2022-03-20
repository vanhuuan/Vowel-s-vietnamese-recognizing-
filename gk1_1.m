close all; clear all;
[data, fs] = audioread('lab_female.wav');
% chuan hoa bien do
data = data / abs(max(data));

% Chia khung
f_d = 0.025;
frames = framing(data, fs, f_d);

% tinh nang luong,ZCR tung khung
[r,c] = size(frames);
ste = 0;
for i = 1 : r
    x = frames(i, :);
    ste(i) = sum(x.^2);  
    ZCR(i) = sum(abs(sign(x(1 : end - 1)) - sign(x(2:end))));
end

% chuan hoa bien do STE
ste = ste./max(ste); 
f_size = round(f_d * fs);
ste_wave = 0;
for j = 1 : length(ste)
    l = length(ste_wave);
    ste_wave(l : l + f_size) = ste(j);
end

%chuan hoa bien do ZCR
ZCRr = ZCR/(2*length(x));
ZCRr = ZCRr/max(ZCRr); 
f_size = round(f_d * fs);
zcr_wave = 0;
for j = 1 : length(ZCRr)
    l = length(zcr_wave);
    zcr_wave(l : l + f_size) = ZCRr(j);
end

%phan doan
T_ZCR=mean(zcr_wave)*1.3;
z=zeros(1,length(data));
for i=1:length(ste_wave)
    if (ste_wave(i)>0.01)
        z(i)=1;
    elseif (zcr_wave(i)<T_ZCR && ste_wave(i)>0.008)
        z(i)=1;
    end
end

%Ve do thi
t = 0 : 1/fs : length(data)/fs; % time (s)
t = t(1:end - 1);
plot(t,data'); xlabel("Time(s)"); hold on;

% ve do thi STE va signal
t1 = 0 : 1/fs : length(ste_wave)/fs;
t1 = t1(1:end - 1);
plot(t1,ste_wave,'r','LineWidth',1.5);

% ve do thi ZCR tren signal
t2 = 0 : 1/fs : length(zcr_wave)/fs;
t2 = t2(1:end - 1);
plot(t2,zcr_wave,'b','LineWidth',1.5); 

%Ve duong bien
for i=1:(length(data)-1)
    if z(i)==0&&z(i+1)==1
        line([i/fs i/fs], [-1 1], 'Color', '#A2142F','LineStyle','--');
    elseif z(i)==1&&z(i+1)==0
        line([i/fs i/fs], [-1 1], 'Color', '#EDB120','LineStyle','--');
    end
end

legend('Speech Signal','STE', 'ZCR');

function [frames] = framing(x,fs,f_d)

% x: signal
% fs: Tan so mau
% f_d: thoi luong khung (in sec)
% frames: khung

f_size = round(f_d * fs);  % Kich thuoc khung
l_s = length(x);    % do dai tin hieu
n_f = floor(l_s/f_size); % so luong khung

% Tao khung
temp = 0;
for i = 1 : n_f
    frames(i,:) = x(temp + 1 : temp + f_size);
    temp = temp + f_size;
end

end

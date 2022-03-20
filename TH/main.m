clear;
close all;

%Cal("01MDA");  
%Cal("07FTC");
%Cal("12FTD");
Cal("16FTH");
function Cal(input_audio)
    nfft = 2048; % số điểm lấy FFT
    draw(input_audio, nfft);
end
function [] = draw(name, nfft)
    % Ham ve do thi nhan doi do: Tin hieu ra, tin hieu vao, nguong loc
    % Ham ste, Matran silence,so frame du, do dai moi frame, 
    % so mau, so frame, Fs
    wbwl = 0.005; % wideband window length
    wbws = 0.003; % wideband window step
    figure('Name',name);
    subplot(5,1,1);
    hold on;
    file = "..\NguyenAmHuanLuyen-16k\"+name+"\a.wav";
    [y,Fs] = audioread(file);
    spectrogram(y, floor(wbwl*Fs), floor((wbwl-wbws)*Fs), nfft, Fs, 'yaxis');
    ylim([0 8]);
    title("Wideband spectrogram of /a/");
    hold off;
    subplot(5,1,2);
    hold on;
    file = "..\NguyenAmHuanLuyen-16k\"+name+"\e.wav";
    [y,Fs] = audioread(file);
    spectrogram(y, floor(wbwl*Fs), floor((wbwl-wbws)*Fs), nfft, Fs, 'yaxis');
    ylim([0 8]);
    title("Wideband spectrogram of /e/");
    hold off;
    subplot(5,1,3);
    hold on;
    file = "..\NguyenAmHuanLuyen-16k\"+name+"\i.wav";
    [y,Fs] = audioread(file);
    spectrogram(y, floor(wbwl*Fs), floor((wbwl-wbws)*Fs), nfft, Fs, 'yaxis');
    ylim([0 8]);
    title("Wideband spectrogram of /i/");
    hold off;
    subplot(5,1,4);
    hold on;
    file = "..\NguyenAmHuanLuyen-16k\"+name+"\o.wav";
    [y,Fs] = audioread(file);
    spectrogram(y, floor(wbwl*Fs), floor((wbwl-wbws)*Fs), nfft, Fs, 'yaxis');
    ylim([0 8]);
    title("Wideband spectrogram of /0/");
    hold off;
    subplot(5,1,5);
    hold on;
    file = "..\NguyenAmHuanLuyen-16k\"+name+"\u.wav";
    [y,Fs] = audioread(file);
    spectrogram(y, floor(wbwl*Fs), floor((wbwl-wbws)*Fs), nfft, Fs, 'yaxis');
    ylim([0 8]);
    title("Wideband spectrogram of /u/");
    hold off;
end
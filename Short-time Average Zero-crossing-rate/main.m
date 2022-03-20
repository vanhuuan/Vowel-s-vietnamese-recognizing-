clear all;
%% Input Speech
[data, fs] = audioread('studio_male.wav');
%[data, fs] = wavread('studio_female.wav');
%[data, fs] = wavread('lab_male.wav');
%[data, fs] = wavread('lab_female.wav');
% Chuan hoa du lieu
data = normalized_data(data);
%% Divide Frame
frame_duration = 0.02;
frame_length = frame_duration * fs;
%chia frame
frames = framing_function(data, fs, frame_length);

%% Zero Crossing Rate
%tan suat di qua 0 cua moi frame
ZCR = zero_crossing_rate_function(frames);
%Chuan hoa du lieu
ZCR = normalized_data(ZCR);

%% Short-Time Energy
%nang luong moi frame
STE = short_time_energy_function(frames);
%Chuan hoa du lieu
STE = normalized_data(STE);
%Gop cac frame lai de ve bieu do
zcr_sample = combine_frame_func(ZCR, frame_length);
ste_sample = combine_frame_func(STE, frame_length);

%% Plot the data,ZCR with STE

subplot(2,1,1);
t = [0 : 1/fs : length(data)/fs];
t = t(1:end - 1);
t1 = [0 : 1/fs : length(zcr_sample)/fs];
t1 = t1(1:end - 1);
t2 = [0 : 1/fs : length(ste_sample)/fs];
t2 = t2(1:end - 1);
plot(t,data,'Color','green');hold on;grid on;
plot(t1,zcr_sample,'r');
plot(t2,ste_sample,'b');
plot(xlim, [1 1]*0.6, '--r');
plot(xlim, [1 1]*0.01, '--b');
hold off;
legend({'Tin hieu vao', 'Zero Crossing Rate', 'Short time energy', 'Nguong ZCRct', 'Nguong Ect'}, 'Location','southeast');
title('Short-time-energy va Zero-crossing-rate');
xlabel('Thoi gian(s)');
%% Xac dinh vi tri tieng noi va khoang lang
[row,col] = size(frames);
%dat nguong
ZCRct = 0.8; Ect = 0.01;
index_zcr = find(ZCR < ZCRct);%loc ra index cua cac mau < nguong ZCRct
index_ste = find(STE > Ect);%loc ra index cua cac mau > nguong Ect
n_frame = 2 * floor(length(data)/frame_length) - 1;
index = [];
for i = 1:n_frame %chay qua cac frame
    j = 1;
    while j <= length(index_zcr)%chay qua index_zcr
        if(i == index_zcr(j))
            k = 1;
            while k <= length(index_ste)%chay qua index_ste
                if(i == index_ste(k))
                    index = [index i];%loc ra index cac frame thoa man < ZCRct va > Ect
                end
                k=k+1;
            end
        end
        j=j+1;
    end
end

%% Plot the Voice and Unvoice
subplot(2,1,2)
times = [0 : 1/fs : col/fs];
times = times(1:end - 1);
plot(t,data,'Color', 'green');hold on;grid on;
% Tim bien thoi gian giao nhau giua tieng noi va khoang lang
index_speech=index(1)-1;
n=2;
for i = 2:length(index)
    timeStart = frame_duration*(index(i)+1)/2;
    timeEnd = frame_duration*(index(i-1)+1)/2;
    if(timeStart - timeEnd >= 0.2)
        index_speech(n)=index(i-1);
        index_speech(n+1)=index(i)-1;
        n=n+2;
    end
end
index_speech(n)=index(i);
local_speech=frame_duration*(index_speech + 1)/2;
% ve do thi phan doan tieng noi va khoang lang
y = [-0.5: 0.5];
for i = 1:length(local_speech)
    if(rem(i,2) == 0)
        plot(local_speech(i)*ones(size(y)), y,'r', 'LineWidth', 1);
    else
        plot(local_speech(i)*ones(size(y)), y,'b', 'LineWidth', 1);
    end
end
hold off
legend('Tin hieu ban dau','Vi tri bat dau co tieng noi','Vi tri ket thuc tieng noi','Location','southeast');
title('Phan doan tieng noi va khoang lang');
xlabel('Thoi gian(s)')

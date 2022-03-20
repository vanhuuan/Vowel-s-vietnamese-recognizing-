clear;
close all;

%Cal("01MDA.wav");  
%Cal("02FVA.wav");
%Cal("03MAB.wav");
%Cal("06FTB.wav");
%Cal("45MDV.wav");
%Cal("44MTT.wav");
%Cal("42FQT.wav");
Cal("30FTN.wav");

function Cal(input_audio)
    [y,t,Fs] = audioread2(input_audio);
    nfft = 512; % số điểm lấy FFT
    %draw(y, t, Fs, input_audio, nfft);
    %Chia khung
    [Frames,Frame_length,NumOfFrame] = seftFraming(y,Fs);
    %Phan doan speech silence
    [Silence,ma,threshold] = Discriminate(Frames,Frame_length,NumOfFrame);
    formants =  findFormant(y, Fs, nfft, Frames, Silence);
    
end
function [y,t,Fs] = audioread2(file)
    % Hàm đọc tín hiệu từ file 
    % Hàm trả về 3 tham số y, t, Fs là tín hiệu đã chuẩn hoá, vector thời
    % gian của tín hiệu và tần số của tín hiệu
    [y,Fs] = audioread(file); % Đọc tín hiệu vào từ file 
    max_value = max(abs(y));
    y = y/max_value; % Chuẩn hoá tín hiệu về từ 0-1 
    t = 1/Fs:1/Fs:(length(y)/Fs); % Vector thời gian của tín hiệu
end
function [Silence,ma,threshold] = Discriminate(Frames,Frame_length,NumOfFrame)
    % Hàm phân biệt khoảng lặng và nguyên âm
    % Hàm nhận vào các Frames của tín hiệu và số Frames được chia
    % Hàm trả về ngưỡng threshold để phân biệt , ma tính được và khoảng
    % lặng Silence tìm được
    [ma] = getMA(Frames,NumOfFrame,Frame_length); % Tính MA của từng khung tín hiệu
    [threshold] = findThreshold(ma); % Tìm threshold từ MA
    dif = @(x) (threshold - ma(x)); % tao han dif voi gia tri duong la khoang lang, am la tieng noi    
    for i = 2:NumOfFrame
        zx(i) = dif(i)*dif(i-1); 
    end  
    Silence = find(zx <= 0); % neu tich 2 phan tu lien tiep mang gia tri am thi o do chinh la bien do can tim, itsect chua gia tri la vi tri cac bien do
end
function [ma] = getMA(frames, n_frames, frame_length)
    % Hàm tính MA của các Frames
    % Hàm trả về MA tính được
    ma = []; % Khởi tạo mảng MA
    for k = 1:n_frames
        sum = 0;
        for j=1:frame_length-1
            sum = sum + abs(frames(j,k)); %tổng các phần tử |x[n-m]|
        end
        ma(k) = sum; %MA[k] = tổng các phần tử
    end
    ma_max = max(ma); % tìm giá tị lớn nhất của MA
    ma = ma.*(1/ma_max); % chuẩn hoá MA
end
function [threshold] = findThreshold(ma) 
    % Hàm tìm ngưỡng Threshold của giản đồ năng lượng MA
    % Hàm nhận vào giản đồ năng lượng MA
    % Hàm trả về threshold
    f = ma; % vecto f để tìm khoảng lặng
    g = ma; % vecto g để tìm nguyên âm
    Tmin = min(ma);
    Tmax = max(ma);
    T = (Tmin+Tmax)/2;
    i = length(f<T); % số lượng giá trị nhỏ hơn T
    p = length(g>T); % số lượng giá trị lớn hơn T
    j = -1;
    q = -1;
    while( j~=i || q~=p) % số lượng không đổi thì dừng vòng lặp
         sumf = 0; % tổng các phần từ nằm trong khoảng lặng
         sumg = 0; % tổng các phần từ nằm trong nguyên âm
         for k = 1:length(f)
             sumf = sumf + max(f(k) - T,0); %tính tổng
             sumg = sumg + max(T - g(k),0); %tính tổng
         end
         if(sumf > sumg) %nếu tổng các phần từ khoảng lặng > tổng các phần tử tiếng nói 
             Tmin = T;        % thì T min = T
         else
             Tmax = T;
         end
         T = (Tmin+Tmax)/2; % Thay đổi T
         j = i;
         q = p;
         i = length(f<T); % tính lại số lượng giá trị nhỏ hơn T
         p = length(g>T); % tính lại số lượng giá trị lớn hơn T
        end
    threshold = (T + Tmin)/2;
end
function [Frames,Frame_length,NumOfFrame] = seftFraming(y0,Fs)  
    % Hàm tách tín hiệu thành các khung tín hiệu dài 20ms
    % Hàm nhận tham số vào là y0, Fs là biên độ và tần số của tín hiệu
    % Hàm trả về y là ma trận chứa biên độ của các khung tín hiệu trong đó
    % Mỗi cột là 1 khung tín hiệu và NumOfFrame là số khung tín hiệu
    Frame_length = floor(0.02*Fs); % tính số điểm có trong 1 khung tín hiệu 20ms 
    Frame_length_2 = floor(0.01*Fs);
    L = length(y0); % xác định độ dài của y0;
    Frames = []; % khởi tạo ma trận rỗng (0x0)  
    a = 1;
    b = Frame_length;  % xác định vị trí đầu và cuối của khung tín hiệu 20ms đầu tiên
    while b<=L     
        Frames = horzcat(Frames,y0(a:b)); % nối mỗi khung tín hiệu vào cột
        a = b-Frame_length_2;
        b = a+Frame_length-1; % xác định vị trí đầu và cuối của khung tín hiệu 20ms tiếp theo
    end
    NumOfFrame = length(Frames(1,:)); % xác định số khung tín hiệu 
    Frames = reshape(Frames, Frame_length, []);
end
function [] = draw(y , t, Fs, name, nfft)
    % Ham ve do thi nhan doi do: Tin hieu ra, tin hieu vao, nguong loc
    % Ham ste, Matran silence,so frame du, do dai moi frame, 
    % so mau, so frame, Fs
    figure('Name',name);
    subplot(2,1,1);
    hold on;
    title('Tin hieo vao');
    xlabel('Thoi gian (s)');
    ylabel('Do lon (dB)');
    plot(t, y, 'b');
    legend('Tin hieu');
    hold off;
    subplot(2,1,2);
    hold on;
    wbwl = 0.005; % wideband window length
    wbws = 0.003; % wideband window step
    spectrogram(y, floor(wbwl*Fs), floor((wbwl-wbws)*Fs), nfft, Fs, 'yaxis');
    ylim([0 8]);
    title("Wideband spectrogram");
    hold off;
end
function formants = findFormant(y, Fs, nfft, Frames, Silence)
    formants = zeros(3,5);
    for i = 1:5
        center = floor((Silence(i*2)+Silence(i*2+1))/2); % Khung ở giữa đoạn nguyên âm
        frame = Frames(:,center);
        window = hamming(length(frame)); %Hàm Hamming window 
        sig = abs(fft(frame.*window));
        sig = log(sig);
        sig = ifft(sig);
        window = cceps(length(sig)); % Tới đoạn này còn đúng 
        sig = abs(fft(frame.*window));
        peak = findpeaks(sig);
        formants(i,1) = peak(1); %f1,f2,f3
        formants(i,2) = peak(2);
        formants(i,3) = peak(3);
    end
    hz5000=5000*length(sig)/Fs;
    f=(0:hz5000)*Fs/length(sig);
    plot(f,sig(1:length(f)));
end
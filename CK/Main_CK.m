clear;
close all;

% Bộ dữ liệu huấn luyện
t1 = [0.45; 0.81; 1.53; 1.85; 2.69; 2.86; 3.78; 4.15; 4.84; 5.14; 135.5; 5.4];
t2 = [0.83; 1.37; 2.09; 2.60; 3.57; 4; 4.76; 5.33; 6.18; 6.68; 239.7; 5.6];
t3 = [1.03; 1.42; 2.46; 2.8; 4.21; 4.52; 6.81; 7.14; 8.22; 8.5; 115.0; 4.5];
t6 = [1.52; 1.92; 3.91; 4.35; 6.18; 6.6; 8.67; 9.14; 10.94; 11.33; 202.9; 15.5];
%Bộ dữ liệu kiểm thử
k30 = [0.59; 0.97; 1.76; 2.11; 3.44; 3.77; 4.70; 5.13; 5.96; 6.28; 233.2; 11.6];
k42 = [0.46; 0.99; 1.56; 2.13; 2.51; 2.93; 3.79; 4.38; 4.77; 5.22; 242.7; 8.5];
k44 = [0.93; 1.42; 2.59; 3; 4.71; 5.11; 6.26; 6.66; 8.04; 8.39; 125.7; 8.5];
k45 = [0.88; 1.34; 2.35; 2.82; 3.76; 4.13; 5.04; 5.5; 6.41; 6.79; 177.8; 5.7];

%Cal("01MDA.wav",t1);  
%Cal("02FVA.wav",t2);
%Cal("03MAB.wav",t3);
%Cal("06FTB.wav",t6);

Cal("45MDV.wav",k45);
Cal("44MTT.wav",k44);
Cal("42FQT.wav",k42);
Cal("30FTN.wav",k30);

function Cal(input_audio,R)
    [y,t,Fs] = audioread2(input_audio);
    n_samples = length(y); %so luong mau cua vecto y

    Frame_duration = 0.02; % Moi khung 20ms
    centers_length = 0.01; % Moi tam cach nhau 10ms
    %Chia khung
    [Frames,Frame_length,NumOfFrame] = SeftFraming(y,Fs);
    %Phan doan speech silence
    [Silence,ma,threshold] = Discriminate(Frames,Frame_length,NumOfFrame);
    %HPS(Frames(:,55),Fs);
    [F0,F0mean,F0std] = Find_F0(Frames,NumOfFrame,Fs,Silence);
    %Ve do thi
    draw(y,threshold,t,ma,Silence,Frame_duration, n_samples,NumOfFrame,Fs,input_audio,centers_length,F0,F0mean,F0std,R);
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
function [Frames,Frame_length,NumOfFrame] = SeftFraming(y0,Fs)  
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
function [] = draw(y , threshold,t, ma, silence,frame_duration, n_samples, n_frames, Fs,name,centers_length,F0,F0mean,F0std,R)
    % Ham ve do thi nhan doi do: Tin hieu ra, tin hieu vao, nguong loc
    % Ham ste, Matran silence,so frame du, do dai moi frame, 
    % so mau, so frame, Fs
    figure('Name',name);
    subplot(5,1,1);
    hold on;
    title('Su dung MA de tach am thanh va khoang lang');
    xlabel('Khung');
    ylabel('Gia tri nang luong da chuan hoa');

    k = 1:n_frames; % k la so khung
    for i = 1:n_frames
        line(i) = threshold;
    end
    plot(k, line, 'r','Linewidth',1); % Ve nguong di qua gian do MA 
    plot(k, ma, 'g', 'Linewidth',1); % Ve gian do nang luong MA
    % Ve cac duong bien
    for i = 2: length(silence)
        j = silence(i);
        plot([j, j], [0, 1], 'b','LineStyle','-.')
    end 
    legend('Nguong','MA', 'Duong bien');
    hold off;
    
    subplot(5,1,2);
    hold on;
    xlabel('Thoi gian thuc (s)');
    ylabel('Bien do tin hieu');
    n = 1:n_samples;
    plot(t, y,'g');
    for i = 2:length(silence)
        if silence(i) > 0
           j = silence(i)*frame_duration/2-centers_length;
           plot([j, j], [-1, 1], 'r','LineStyle','-.')
        end
    end
    for i = 1:length(R)-2
        j = R(i);
        plot([j, j], [-1, 1], 'b','LineStyle','--');
    end
    legend('Tin hieu','Bien tu tim kiem','Bien Chuan');
    hold off;
    subplot(5,1,3);
    hold on;
    plot(k, F0,'.');
    for i = 2:length(silence)
        j = silence(i);
        plot([j, j], [0,300], 'r','LineStyle','--');
    end
    title("Duong bieu dien tan so co ban" + ' : F0mean= '+F0mean + ' F0std= '+F0std ); 
    xlabel("Khung"); ylabel("F_0 (Hz)");
    legend('F0','Đường biên');
    hold off;
    stdChuan = R(length(R));
    meanChuan = R(length(R)-1);
    saisoStd = ((F0std-stdChuan)/stdChuan)*100
    saisoMean = ((F0mean - meanChuan)/meanChuan) * 100
    subplot(5,1,4);
    hold on;
    f= y(0.02*Fs:0.04*Fs);
    nfft = 2048;
    window = hamming(length(f)); %Hamming window function
    FFT = abs(fft(f.*window,nfft));
    FFT = FFT(1:length(FFT)/2);
    plot(FFT)
    title("Frame không tuần hoàn"); 
    hold off;
    subplot(5,1,5);
    hold on;
    f= y(0.90*Fs:0.92*Fs);
    FFT = abs(fft(f.*window,nfft));
    FFT = FFT(1:length(FFT)/2);
    plot(FFT)
    title("Frame tuần hoàn"); 
    hold off;
end
function [F0] = CalculateF0(Frame,Fs)
    % Hàm tính của một Frame
    % Hàm nhận vào 1 Frame và tần số lấy mẫu của tín hiệu
    window = hamming(length(Frame)); %Hàm Hamming window 
    n = 2048; % số điểm lấy FFT
    FFT = abs(fft(Frame.*window,n));
    FFT = FFT(1:length(FFT)/2); % 0 -> pi, pi -> 2pi là như nhau nên cắt bỏ 1 nửa
    l = Fs/n; % Độ phân giải của phổ
    distance_peaks = 4; %Khoảng cách giữa các peak
    [y, x] = findpeaks(FFT,'MinPeakDistance',distance_peaks);
    num = 5; % Số đỉnh cần tìm
    distance = [];
    for i = 1:num-1
        distance(i) = x(i+1) - x(i);
    end
    F0 = max(distance)*l; %F0 bằng khoảng cách giữa 2 điểm cạnh nhau có độ dài lớn nhất
end
function [F0,F0mean,F0std] = Find_F0(Frames,n_frames,Fs,Silence)
    % Hàm tìm F0 của các Frame đã chia nguyên âm và khoảng lặng.
    % Hàm nhận vào các mốc của khoảng lặng Silence, các Frames , số lượng
    % frame n_frame và tần số lấy mẫu Fs.
    % Hàm trả về đường F0 , F0 mean và F0 std.
    F0 = []; % Ma trận chứa f0
    raw_f0 = []; % Ma trận chứa f0 chưa qua xử lí bao gồm tuần hoàn và không tuần hoàn
    for i = 1:n_frames
        raw_f0 =  [raw_f0 CalculateF0(Frames(:,i),Fs)]; % Tính f0 của từng frame
    end
    si = false;
    j=1;
    for i = 1:n_frames % chia các frame trong khoảng lặng có f0 = 0
        if  j<= length(Silence) && i == Silence(j)
            j=j+1;
            si = not(si);
        end
        if not(si)
            F0(i) = raw_f0(i);
        else
            F0(i) = 0;
        end
    end
    [F0std,F0mean] = FindF0mean_F0std(F0);
    canTren = F0std + F0mean;
    canDuoi = F0mean-F0std;
    for i = 1:n_frames % Xét tuần hoàn hay không tuần hoàn
        if F0(i) >= canDuoi && F0(i) <= canTren && F0(i) >= 70 && F0(i) <= 400 %F0 nằm trong khoảng từ cận duưới đến cận trên của tín hiệu và nằm trong khoản tai người nghe đc
            F0(i) = raw_f0(i);
        else
            F0(i) = 0;
        end
    end
    [F0std,F0mean] = FindF0mean_F0std(F0) %Tính F0mean va F0 std
end
function [F0std,F0mean] = FindF0mean_F0std(F0) 
    % Hàm tìm giá trị F0 mean và F0 std
    % Hàm nhận tham số đầu vào F0 là đường F) đã tính
    % Hàm trả về giá trị F0 mean và F0 std
    F0_new = []; % Lưu giá trị khác 0 của F0
    l = 1:length(F0); % Xác định độ dài F0
    for i=l % Lọc các giá trị bằng 0
        if (F0(i)~=0)         
            F0_new = [F0_new F0(i)]; % Nếu F0 khác không thì thêm vào F0_new 
        end
    end
    F0mean = mean(F0_new);
    F0std = std(F0_new);
end
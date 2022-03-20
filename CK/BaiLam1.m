clear;5.
close all;

t1 = [0.45; 0.81; 1.53; 1.85; 2.69; 2.86; 3.78; 4.15; 4.84; 5.14; 135.5; 5.4];
t2 = [0.83; 1.37; 2.09; 2.60; 3.57; 4; 4.76; 5.33; 6.18; 6.68; 239.7; 5.6];
t3 = [1.03; 1.42; 2.46; 2.8; 4.21; 4.52; 6.81; 7.14; 8.22; 8.5; 115.0; 4.5];
t6 = [1.52; 1.92; 3.91; 4.35; 6.18; 6.6; 8.67; 9.14; 10.94; 11.33; 202.9; 15.5];

k30 = [0.59; 0.97; 1.76; 2.11; 3.44; 3.77; 4.70; 5.13; 5.96; 6.28; 233.2; 11.6];
k42 = [0.46; 0.99; 1.56; 2.13; 2.51; 2.93; 3.79; 4.38; 4.77; 5.22; 242.7; 8.5];
k44 = [0.93; 1.42; 2.59; 3; 4.71; 5.11; 6.26; 6.66; 8.04; 8.39; 125.7; 8.5];
k45 = [0.88; 1.34; 2.35; 2.82; 3.76; 4.13; 5.04; 5.5; 6.41; 6.79; 177.8; 5.7];

%Cal("01MDA.wav",t1);  
%Cal("02FVA.wav",t2);
%Cal("03MAB.wav",t3);
%Cal("06FTB.wav",t6);

%Cal("45MDV.wav",k45);
%Cal("44MTT.wav",k44);
%Cal("42FQT.wav",k42);
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
    draw(y,y,Frames,threshold,t,ma,Silence,Frame_duration, n_samples,NumOfFrame,Fs,input_audio,centers_length,F0,F0mean,F0std,R);
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
    [ma] = getMA(Frames,NumOfFrame,Frame_length);
    [threshold] = findThreshold(ma);
    dif = @(x) (threshold - ma(x)); % tao han dif voi gia tri duong la khoang lang, am la tieng noi    
    for i = 2:NumOfFrame
        zx(i) = dif(i)*dif(i-1); 
    end  
    Silence = find(zx <= 0); % neu tich 2 phan tu lien tiep mang gia tri am thi o do chinh la bien do can tim, itsect chua gia tri la vi tri cac bien do
end
function [ma] = getMA(frames, n_frames, frame_length)
    % Ham nhan cac doi sola : frames : vector frames,n_frames: do dai cua vector am thanh
    % goc, frame_length: so mau cua 1 frame
    % chieu dai khung
    % Tra ve gian do nang luong MA
    ma = [];
    for k = 1:n_frames
        sum = 0;
        for j=1:frame_length-1
            sum = sum + abs(frames(j,k)); %tong cac phan tu |x[n-m]|
        end
        ma(k) = sum; %MA[k] = tong cac phan tu
    end
    ma_max = max(ma); %tim gia tri lon nhat cua MA
    ma = ma.*(1/ma_max); %chuyen doi bien do cua MA lon nhat thanh 1.
end % Ket thuc ham getMA
function [threshold] = findThreshold(ma) % ham tim threshold cua gian do nang luong MA
    %[a,idx]=unique(ma); % tim overlap cua gian do nang luong MA
    %[ii,ii]=sort(idx); % a(ii) la overlap cua gian do nang luong MA
    f = ma; % vecto f de tim khoang lang
    g = ma; % vecto g de tim giong noi
    Tmin = min(ma);
    Tmax = max(ma);
    T = (Tmin+Tmax)/2;
    i = length(f<T); % so luong gia tri nho hon T
    p = length(g>T); % so luong gia tri lon hon T
    j = -1;
    q = -1;
    while( j~=i || q~=p) % so luong khong doi thi dung vong lap
         sumf = 0; % tong cac phan tu nam trong khoang lang
         sumg = 0; % tong cac phan tu nam trong tieng noi
         for k = 1:length(f)
             sumf = sumf + max(f(k) - T,0); %tinh tong
             sumg = sumg + max(T - g(k),0); %tinh tong
         end
         if(sumf > sumg) %neu tong cac phan tu khoang lang > tong cac phan tu tieng noi 
             Tmin = T;        % thi T min = T
         else
             Tmax = T;
         end
         T = (Tmin+Tmax)/2; % Thay doi T
         j = i;
         q = p;
         i = length(f<T); % tinh lai so luong gia tri nho hon T
         p = length(g>T); % tinh lai so luong gia tri lon hon T
        end
    threshold = (T + Tmin)/2;
end
function [] = draw(processedSignal,y , frames, threshold,t, ma, silence,frame_duration, n_samples, n_frames, Fs,name,centers_length,F0,F0mean,F0std,R)
    % Ham ve do thi nhan doi do: Tin hieu ra, tin hieu vao, nguong loc
    % Ham ste, Matran silence,so frame du, do dai moi frame, 
    % so mau, so frame, Fs
    figure('Name',name);
    subplot(3,1,1);
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
    
    subplot(3,1,2);
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
    subplot(3,1,3);
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
end
function [F0_HPS] = HPS(Frame,Fs)
    window = hamming(length(Frame)); %Hamming window function
    n = 2048; % So diem lay FFT
    FFT = abs(fft(Frame.*window,n));
    FFT = FFT(1:length(FFT)/2); % 0 -> pi, pi -> 2pi là nh
    %figure(1);
    %plot(FFT);
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
    [x,y]=findpeaks(hps, 'SORTSTR', 'descend');
    Maximum = y(1);
    F0_index = 1;
    F0_next_index=2;
    while true
        Maximum = y(F0_index);
        if y(F0_next_index) > Maximum*0.29 && y(F0_next_index) < Maximum*0.60
            F0_index = F0_next_index;
            F0_next_index = F0_next_index + 1;
            continue;
        else
            break;
        end
    end
    F0_HPS = (Maximum/n)*Fs;
end
function [F0,F0mean,F0std] = Find_F0(Frames,n_frames,Fs,Silence)
    F0 = []; %Ma tran chua f0
    hps_f0 = [];
    for i = 1:n_frames
        hps_f0 =  [hps_f0 HPS(Frames(:,i),Fs)];
    end
    si = false;
    j=1;
    for i = 1:n_frames
        if  j<= length(Silence) && i == Silence(j)
            j=j+1;
            si = not(si);
        end
        if not(si)
            F0(i) = hps_f0(i);
        else
            F0(i) = 0;
        end
    end
    [F0std,F0mean] = FindF0mean_F0std(F0);
    canTren = F0std + F0mean;
    canDuoi = F0std - F0mean;
    for i = 1:n_frames
        if F0(i) >= canDuoi && F0(i) <= canTren && F0(i) >= 70 && F0(i) <= 400 %F0 nằm trong khoảng từ cận duưới đến cận trên của tín hiệu và nằm trong khoản tai người nghe đc
            F0(i) = hps_f0(i);
        else
            F0(i) = 0;
        end
    end
    [F0std,F0mean] = FindF0mean_F0std(F0);
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
function y = MedSmoothing(x,N)
    % Hàm lọc trung vị
    % Hàm trả về y là tín hiệu sau khi lọc
    % Hàm nhận x là tín hiệu trước khi lọc, N la bậc cần lọc
    y = x; %khoi tao gia tri tin hieu y
    temp = 1:N; % tao truc thoi gian tin hieu cua so co chieu dai N
    for k = 1:length(x) % k chay tu 1 den chieu dai cua tin hieu x
        for i = 1:N % i chạy từ 1 đến N
            if(k < ceil(N/2))
                % xét giá trị tín hiệu của temp bên trái tín hiệu x
                if(i < ceil(N/2)-k+1) % i < phần tử giữa của temp
                    temp(i) = 0;
                else
                    temp(i) = x(k+i-ceil(N/2)); % i > phần tử giữa của temp
                end
            else
                % xét giá trị tín hiệu của temp bên phải tín hiệu x
                if(k > length(x) - ceil(N/2) + 1)
                    %i < phần tử giữa của temp 
                    if(i < length(x) + ceil(N/2) - k + 1)
                        temp(i) = x(k-ceil(N/2) + i);
                    else
                        temp(i) = 0; % i > phần tử giữa của temp
                    end
                else
                    % Xét tín hiệu temp chạy ở giữa
                    temp(i) = x(k-ceil(N/2)+i);
                end
            end
        end
    % Hàm median tìm phần tử trung vị của dãy
    y(k) = median(temp);
    end
end
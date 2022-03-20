clear;
close all;
Cal("phone_F2.wav");  
Cal("phone_M2.wav");
Cal("studio_F2.wav");
Cal("Studio_M2.wav");  

function Cal(file)
    % Hàm dùng để thực hiện từng bước việc tính F0 và hiện F0 ra màn hình.
    % Hàm nhận tham số file là đường dẫn đến file tín hiệu
    % Hàm không trả về giá trị
    figure('Name',file);
    [y,t,Fs] = audioread2(file);% Đọc tín hiệu đầu vào
    [yuv,lag] = AutocorrelationFunction(y(1.54*Fs:1.56*Fs)); % Hàm tự tương quan với đoạn UnVoice
    if file == "phone_F2.wav"
        [yuv,lag] = AutocorrelationFunction(y(1.90*Fs:1.92*Fs)); % Hàm tự tương quan với đoạn UnVoice
    end
    if file == "studio_F2.wav"
        [yuv,lag] = AutocorrelationFunction(y(1.79*Fs:1.81*Fs)); % Hàm tự tương quan với đoạn UnVoice
    end
    subplot(5,1,4); plot(lag,yuv); title("UnVoice"); xlabel("Lag"); % Hiển thị đoạn UnVoice
    [y1,t1] = CutUseLessSignal(y,t,0.1); % tách vùng tín hiệu theo biên độ
    [yv,lag] = AutocorrelationFunction(y(1.02*Fs:1.04*Fs)); % Hàm tự tương quan với đoạn Voice
    subplot(5,1,5); plot(lag,yv); title("Voice"); xlabel("Lag"); % Hiển thị đoạn Voice
    subplot(5,1,1); plot(t,y); hold on; plot(t1,y1); title("Tin hieu ban dau"); xlabel("Thoi gian (s)"); ylabel("Bien do"); 
    subplot(5,1,2); plot(t1,y1); title("Tin hieu sau khi tach theo bien do"); xlabel("Thoi gian (s)");ylabel("Bien do"); 
    [y2,NumOfFrame] = SeftFraming(y1, Fs); % Chia tín hiệu thành các khung có độ dài 20ms 
    [F0,tF0] = FindF0(y2,t1,NumOfFrame,Fs);% Xác định đường F0
    F0 = MedSmoothing(F0,3);% Kết hợp lọc trung vị để loại bỏ pitch ảo  
    [F0std, F0mean] = FindF0mean_F0std(F0); % Xác định F0 mean và F0 std
    subplot(5,1,3); plot(tF0,F0,'.'); title("Duong bieu dien tan so co ban" + ' : F0mean= '+F0mean + ' F0std= '+F0std ); 
    xlabel("Thoi gian (s)"); ylabel("F_0 (Hz)");
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

function [y,NumOfFrame] = SeftFraming(y0,Fs)  
    % Hàm tách tín hiệu thành các khung tín hiệu dài 20ms
    % Hàm nhận tham số vào là y0, Fs là biên độ và tần số của tín hiệu
    % Hàm trả về y là ma trận chứa biên độ của các khung tín hiệu trong đó
    % mỗi cột là 1 khung tín hiệu và NumOfFrame là số khung tín hiệu
    SoDiemTrong1Khung = floor(0.02*Fs); % tính số điểm có trong 1 khung tín hiệu 20ms 
    L = length(y0); % xác định độ dài của y0;
    y = []; % khởi tạo ma trận rỗng (0x0)  
    a = 1;
    b = SoDiemTrong1Khung;  % xác định vị trí đầu và cuối của khung tín hiệu 20ms đầu tiên
    while b<=L     
        y = horzcat(y,y0(a:b)); % nối mỗi khung tín hiệu vào cột
        a = b;
        b = a+SoDiemTrong1Khung-1; % xác định vị trí đầu và cuối của khung tín hiệu 20ms tiếp theo
    end
    NumOfFrame = length(y(1,:)); % xác định số khung tín hiệu 
end

function [Cut_Y, Cut_T] = CutUseLessSignal(y,t,A) 
    % Hàm loại bỏ 2 đoạn silence ở 2 đầu tín hiệu 
    % Hàm nhận vào y,t là biên bộ và thời gian của tín hiệu gốc, A là mức
    % phân biệt 2 đoạn silence
    % Hàm trả về trục biên độ và trục thời gian của tín hiệu vừa cắt
    % (Cut_Y, Cut_T)
    l = length(y);% xac dinh so diem trong y  
    last = l; 
    % Xác định vị trí đầu của đoạn cần tính 
    for i = 1:last     
        if abs(y(i))>A         
            first = i;          
            break;     
        end
    end 
    % Xác định vị trí cuối của đoạn cần tính 
    for i = l:-1:1     
        if abs(y(i))>A         
            last = i;         
            break;     
        end
    end
    Cut_Y = y(first:last);
    Cut_T = t(first:last);
end

function [F0,tF0] = FindF0(y2,t2,NumOfFrame,Fs) 
    % Hàm tìm F0 sử dụng hàm tự tương quan
    % Hàm nhận giá trị đầu vào là biên độ của các khung tín hiệu y2 , số
    % khung NumOfFrame và tần số của tín hiệu Fs
    % Hàm trả về đường F0 và thời gian tương ứng tF0
    F0 = []; % khởi tạo ma trận rỗng (0x0)  
    tF0 = [];
    Peak = [];
    Lag =[];
    for i = 1:NumOfFrame %tìm peak và lag của từng khung tín hiệu
        [y,l] = AutocorrelationFunction (y2(:,i)); % dùng hàm tự tương quan
        [p,l] = Find_Peak(y); % tìm peak và lag của khung
        Peak = [Peak p];
        Lag = [Lag l];
    end
    ThreshHold = mean(Peak)-std(Peak); % ngưỡng của tín hiệu
    t = t2(1);
    min_ = Fs/400; max_ = Fs/70; % giới hạn của lag
    for i = 1:NumOfFrame  % Duyệt từng khung để tìm F0
        F00 = 0; % F0 mặc định nếu dưới ngưỡng
        if Peak(i) > ThreshHold &&  min_<=Lag(i) && Lag(i)<=max_ 
            F00 = Fs/Lag(i); % nếu đạt ngưỡng thì tìm giá trị F0
        end
        F0 = [F0 F00];  % Thêm F0 và thời gian tương ứng
        tF0 = [tF0 t];
        t = t + 0.02;
    end
end

function [y,lag] = AutocorrelationFunction (signal) 
    % Hàm tự tương quan
    % Hàm nhận 1 khing tín hiệu signal
    % Hàm trả về y chứa giá trị hàm tự tương quan và lag tương tứng
    L = length(signal); % xác định độ dài tín hiệu
    lag = 1:L; % xác định trục độ trễ của tín hiệu 
    y = []; % khởi tạo ma trận chứa giá trị của hàm tự tương quan
    for k = lag     
        sum = 0;% % đặt tổng bằng 0    
        for j = 1:L-k         
            B = 0;       
            if j+k<=L             
                B = signal(j+k);         
            end
            sum = sum + signal(j)*B; % áp dụng công thức 
        end
        y = [y sum]; % thêm giá trị tính được vào y
    end
    y = y/max(abs(y));
end

function [peak,lag] = Find_Peak(acf)
    % Hàm dùng thuật toán Maximum autocorrelation peak detection để tìm
    % Peak max và lag tương ứng của hàm tự tương quan
    % Hàm nhận acf là giá trị của hàm tự tương quan
    % Hàm trả về max peak và lag của nó
    p = zeros(1,length(acf)); % Lưu peak point
    l = zeros(1,length(acf)); % Lưu lag
    for i = 2:1:(length(acf)-1) % Tìm max Peak và lag của nó
        if(acf(i)>acf(i-1) && acf(i)>acf(i+1))
           p(i) = acf(i);
           l(i) = i;
        end  
    end 
    peak = -1; % max peak
    lag = -1; % lag tương ứng
    for i = 2:1:(length(acf)-1) % Tìm max peak và lag của nó
        if p(i) > peak
            peak = p(i);
            lag = l(i);
        end
    end
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
close();
%[y,Fs] = Cal("phone_F1.wav");
[y,Fs] = audioread("phone_F1.wav"); 
[F0] = CalF0Use_AutocorrelationFunction (y(1.37*Fs:1.39*Fs),Fs)
[F0] = CalF0Use_AutocorrelationFunction (y(1.00*Fs:1.02*Fs),Fs)
%plot(lag,y);

function [F0,t1] = Cal(file)
    [y,t,Fs] = audioread2(file);% doc tin hieu dau vao 
    [y1,t1] = CutUseLessSignal(y,t,0.1); %tach vung tin hieu theo bien do 
    [y2,NumOfFrame] = SeftFraming(y1, Fs); %tach vung tinh hieu thanh nhieu khung tin hieu 30ms 
    [F0,tF0] = FindF0(y2,t1,NumOfFrame,Fs);%xac dinh duong F0 
    F0 = MedSmoothing(F0,3);%ket hop loc trung vi de lam tron duong F0  
    [F0std, F0mean] = FindF0mean_F0std(F0); % xac dinh RMSE va tan so co ban F0 trung binh
end
function [y,t,Fs] = audioread2(file) % doc tin hieu doc vao tu file tra ve 3 tham so y, t, Fs 
    [y,Fs] = audioread(file); % doc tin hieu doc vao tu file 
    max_value = max(abs(y)); % tim bien do lon nhat cua tin hieu y 
    y = y/max_value; % quy bien do cua tin hieu cao nhat ve 1 
    t = 1/Fs:1/Fs:(length(y)/Fs); % thoi gian cua y
end

function [y,NumOfFrame] = SeftFraming(y0,Fs)  
    SoDiemTrong1Khung = floor(0.02*Fs); % so diem trong 1 khung tin hieu 30ms 
    L = length(y0); % xac dinh do dai cua y0;
    y = []; % khoi tao ma tran rong (0x0)  
    a = 1;
    b = SoDiemTrong1Khung;  % xac dinh vi tri dau va cuoi cua khung tin hieu 30ms dau tien 
    while b<=L     
        y = horzcat(y,y0(a:b)); % noi theo so cot (tuc la tang so cot)        
        a = b;
        b = a+SoDiemTrong1Khung-1; % xac dinh vi tri dau va cuoi cua khung tin hieu 30ms dau tien 
    end
    NumOfFrame = length(y(1,:)); % xac dinh so khung tin hieu  
end

function [Cut_Y, Cut_T] = CutUseLessSignal(y,t,A) 
    l = length(y);% xac dinh so diem trong y  
    last = l; 
    % xac dinh vi tri dau cua khoang can tinh 
    for i = 1:last     
        if abs(y(i))>A         
            first = i;          
            break;     
        end
    end 
    % xac dinh vi tri cuoi cua khoang can tinh 
    for i = l:-1:1     
        if abs(y(i))>A         
            last = i;         
            break;     
        end
    end
    Cut_Y = y(first:last); % cat vung tuan hoan gan cho Y_TuanHoan 
    Cut_T = t(first:last); % cat vung tuan hoan gan cho t_TuanHoan 
end

function [F0,tF0] = FindF0(y2,t1,NumOfFrame,Fs) 
    F0 = []; % khoi tao ma tran rong (0x0)  
    tF0 = [];
    t = t1(1);
    for i = 1:NumOfFrame         
        [F00] = CalF0Use_AutocorrelationFunction(y2(:,i),Fs); % tinh tan so co ban cua tung khung tin hieu 30ms 
        F0 = [F0 F00]; % them tan so co ban tinh duoc vao ma tran  
        tF0 = [tF0 t];
        t = t + 0.02;
    end
end

function [y,lag] = AutocorrelationFunction (signal)  
    L = length(signal); % xac dinh do dai cua tinh hieu 
    lag = 1:L; % xac dinh truc do tre (lag) cua y  
    y = []; %khoi tao ma tran rong (0x0) 
    for k = lag     
        sum = 0;% khoi tao lai gia tri 0     
        for j = 1:L-k         
            B = 0;       
            if 1<=j+k && j+k<=L             
                B = signal(j+k);         
            end
            sum = sum + signal(j)*B;     
        end
        y = [y sum]; % them gia tri tinh duoc vao y 
    end
    y = y/max(abs(y));
    p = Find_Peak(y);
    disp(max(p));
    figure();
    plot(lag,y);
end

function p = Find_Peak(acf) % dung Maximum autocorrelation peak detection
    p = zeros(1,length(acf)); % luu peak point
    for i = 2:1:(length(acf)-1)
        if(acf(i)>acf(i-1) && acf(i)>acf(i+1))
           p(i) = acf(i); 
        end  
    end   
end
function [F0] = CalF0Use_AutocorrelationFunction (Signal,Fs)  
    F0 = 0;  
    [Yacf,lag] = AutocorrelationFunction(Signal); % xac dinh ham tu tuong quan 
    if max(abs(Yacf)) < 0.2 % la nguong co phai vung tuan hoan hay khong    
        return  
    end  
    min_ = Fs/400; max_ = Fs/70; % xac dinh khoang do tre 70 - 400
    L = length(lag);   
    Yacfmax = -1;
    for i = 1:L     
        if Yacf(i)>Yacfmax && min_<=lag(i) && lag(i)<=max_      
            Yacfmax = Yacf(i); % tim vi tri co bien do lon nhat cua ham tu tuong quan               
            F0 = Fs/lag(i); % tinh tan so co ban F0 = 1 / T0, nhan voi Fs la de co cung do lon voi Fs     
        end
    end
end
F0chuan = 131.58; 
[y,t,Fs] = audioread2("phone_female.wav");% doc tin hieu dau vao 
[y1,t1] = TachVungTinHieuTheoBienDo(y,t,0.1); % (BƯỚC 1) tach vung tin hieu theo bien do 
[y2,t2,SoKhung] = TachThanhCacKhungTinHieu(y1, t1, Fs); % (BƯỚC 2) tach vung tinh hieu thanh nhieu khung tin hieu 30ms 
[F0,tF0] = XacDinhDuongF0(y2,t2,SoKhung,Fs);% (BƯỚC 3 và BƯỚC 4) xac dinh duong F0 
F0 = medfilt1(F0,7);% (BƯỚC 5) ket hop loc trung vi de lam tron duong F0  
[RMSE, F0tb] = XacDinhRMSE_TanSoF0(F0,F0chuan) % xac dinh RMSE va tan so co ban F0 trung binh   
% ve minh hoa 
subplot(3,1,1); plot(t,y);
hold on; plot(t1,y1); 
title("Tin hieu ban dau"); 
xlabel("Thoi gian (s)");
ylabel("Bien do"); 
subplot(3,1,2); 
plot(t1,y1); 
title("Tin hieu sau khi tach theo bien do"); 
xlabel("Thoi gian (s)");ylabel("Bien do"); 
subplot(3,1,3); plot(tF0,F0,"*"); 
title("Duong bieu dien tan so co ban"); 
xlabel("Thoi gian (s)"); 
ylabel("F_0 (Hz)");

function [y,t,Fs] = audioread2(file) % doc tin hieu doc vao tu file tra ve 3 tham so y, t, Fs 
    %--------------------------------------- 
    % [y,t,Fs] = audioread2(file) 
    % y: tin hieu sau khi doc file 
    % t: truc thoi gian 
    % Fs: tan so lay mau 
    [y,Fs] = audioread(file); % doc tin hieu doc vao tu file 
    max_value = max(abs(y)); % tim bien do lon nhat cua tin hieu y 
    y = y/max_value; % quy bien do cua tin hieu cao nhat ve 1 
    t = 1/Fs:1/Fs:(length(y)/Fs); % thoi gian cua y
end

function [Y_TuanHoan, t_TuanHoan] = TachVungTinHieuTheoBienDo(y,t,A) 
    % tra ve vung tin hieu tuan hoan 
    % ------------------------------ 
    % Y_TuanHoan = tin hieu trong vung tuan hoan 
    % t_TuanHoan = truc thoi gian cua tin Y_TuanHoan 
    % y = tin hieu dau vao 
    % t = truc thoi gian cua y 
    % A = bien do de xac dinh vung tuan hoan
    dodai = length(y);% xac dinh so diem trong y dau = 1;  
    cuoi = dodai; 
    % xac dinh vi tri dau cua vung tuan hoan 
    for i = 1:dodai     
        if abs(y(i))>A         
            dau = i;          
            break;     
        end
    end
    % --------------------------------------- 
    % xac dinh vi tri cuoi cua vung tuan hoan 
    for i = dodai:-1:1     
        if abs(y(i))>A         
            cuoi = i;         
            break;     
        end
    end
    % --------------------------------------- 
    Y_TuanHoan = y(dau:cuoi); % cat vung tuan hoan gan cho Y_TuanHoan 
    t_TuanHoan = t(dau:cuoi); % cat vung tuan hoan gan cho t_TuanHoan 
end

function [y,t,SoKhung] = TachThanhCacKhungTinHieu(y0,t0,Fs) 
    % ham tra ve tat cac khung tin hieu 30ms, voi cac khung dan xen nhau 
    % ------------------------------------------------------------------ 
    % [y,t,SoKhung] = TachThanhCacKhungTinHieu(y0,t0,Fs) 
    % y = ma tran MxN, voi M: so diem tung khung tin hieu 30ms, N: so khung 
    % t = ma tran NxM, tuong ung la truc thoi gian cua tung khung 30ms 
    % y0 = tin hieu dau vao 
    % t0 = truc thoi gian cua y0 
    % Fs = tan so lay mau cua y0 
    SoDiemTrong1Khung = floor(0.03*Fs); % so diem trong 1 khung tin hieu 30ms 
    L = length(y0); % xac dinh do dai cua y0;
    y = []; % khoi tao ma tran rong (0x0) 
    t = []; % khoi tao ma tran rong (0x0) 
    a = 1;
    b = SoDiemTrong1Khung;  % xac dinh vi tri dau va cuoi cua khung tin hieu 30ms dau tien 
    while b<=L     
        y = horzcat(y,y0(a:b)); % noi theo so cot (tuc la tang so cot)    
        t = vertcat(t,t0(a:b)); % noi theo so hang (tuc la tang so hang)    
        a = b-floor(0.029*Fs);
        b = a+SoDiemTrong1Khung-1; % xac dinh vi tri dau va cuoi cua khung tin hieu 30ms dau tien 
    end
    SoKhung = length(y(1,:)); % xac dinh so khung tin hieu  
    disp(SoKhung)
end

function [y,lag] = HamTuTuongQuan(tinhieu) 
    % ham tra ve ham tu tuong quan cua tin hieu dau vao (tinhieu) 
    % ------------------------------------------------------------ 
    % [y,lag] = HamTuTuongQuan(tinhieu) 
    % y = ham tu tuong quan cua tin hieu dau vao
    % lag = truc do tre (lag) cua y  
    % tinhieu = tin hieu dau vao de tim ham tu tuong quan 
    dodai = length(tinhieu); % xac dinh do dai cua tinh hieu 
    lag = 1:dodai; % xac dinh truc do tre (lag) cua y  
    y = []; %khoi tao ma tran rong (0x0) 
    for k = lag     
        tong = 0;% khoi tao lai gia tri 0     
        for j = 1:dodai         
            B = 0; % vi tri khong xac dinh se cho bang 0         
            if 1<=j+k && j+k<=dodai             
                B = tinhieu(j+k);         
            end
            tong = tong + tinhieu(j)*B;     
        end
        y = [y tong]; % them gia tri tinh duoc vao y 
    end
end

function [F0,tF0] = XacDinhDuongF0(y2,t2,SoKhung,Fs) 
    % xac dinh duong F0 tu tat ca cac khung tin hieu 30ms % tra ve sai so toan phuong trung binh RMSE 
    % -------------------------------------------------------------- 
    % F0 = duong tan so co ban F0 
    % tF0 = truc cua F0 tren do thi 
    % RMSE = sai so toan phuong trung binh 
    % y2 = ma tran MxN, chua tat cac cac khung tin hieu 30ms 
    % t2 = ma tran NxM, chua truc thoi gian tung khung tin hieu 30ms 
    % Fs = tan so lay mau % F0chuan = tan so co ban F0 chuan do thu cong % 
    tF0 = []; % khoi tao ma tran rong (0x0) 
    F0 = []; % khoi tao ma tran rong (0x0)  
    for i = 1:SoKhung         
        [F00] = TinhTanSoCoBan_SDHTTQ(y2(:,i),Fs); % tinh tan so co ban cua tung khung tin hieu 30ms         
        F0 = [F0 F00]; % them tan so co ban tinh duoc vao ma tran         
        tF0 = [tF0 t2(i,floor((1+length(t2(i,:)))/2))]; % xac dinh truc F0 tren do thi  
    end
end

function [RMSE,F0tb] = XacDinhRMSE_TanSoF0(F0,F0chuan) 
    % tra ve RMSE cua duong tan so co ban F0 
    % --------------------------------------- 
    % [RMSE,F0tb] = XacDinhRMSE_TanSoF0(F0,F0chuan) 
    % RMSE = sai so toan phuong trung binh 
    % F0tb = Tan so co ban trung binh 
    % F0 = duong tan so co ban 
    % F0chuan = tan so co ban chuan do thu cong 
    L = 1:length(F0); % xac dinh do dai duong F0 
    Sum = 0; % khoi tao gia tri bang 0 
    Dem=0; % bien dem so luong F0 xac dinh 
    F0tb=0; 
    for i=L     
        if (F0(i)>0)         
            Dem = Dem+1;         
            Sum = Sum + (F0chuan-F0(i))^2; % tong cac binh phuong sai
            F0tb = F0tb+F0(i);     
        end
    end
    RMSE = sqrt(Sum/Dem); % tinh sai so toan phuong trung binh RMSE 
    F0tb = F0tb/Dem;% tinh tan so co ban trung binh
end

function [F0,Yttq,lag,lagmax] = TinhTanSoCoBan_SDHTTQ(TinHieu,Fs) 
    % tinh tan so co ban 
    % tra ve ham tu tuong quan voi cac do tre (lag) khac nhau,  
    % tra ve do tre (lag) tai do dinh cao nhat thoa man 80<=F0<=400 
    % ---------------------------------------------------------------- 
    % F0 = tan so co ban ma thuat toan xac dinh duoc 
    % Yttq = ham tu tuong quan cua tin hieu dau vao 
    % lag = truc do tre (lag) cua Yttq 
    % lagmax = do tre (lag) tai do dinh cao nhat thoa man 80<=F0<=400 
    F0 = 0;
    lagmax = 0; % khoi tao cac gia tri ban dau 
    [Yttq,lag] = HamTuTuongQuan(TinHieu); % xac dinh ham tu tuong quan 
    % kiem tra tin hieu vao co tuan hoan khong? 
    if Yttq(int32((1+length(lag))/2)) <10     
        return  
    end
    %----------------------------------------------------------------- 
    % su dung noi suy cubic spline ----------------------------------- 
    %L = max(lag); 
    %lag2 = -L:0.1:L; 
    %Yttq = spline(lag,Yttq,lag2); 
    %lag = lag2; 
    % ---------------------------------------------------------------- 
    % tim kiem do tre co dinh lon nhat (Fs/400 < dotre <Fs/80) 
    a = Fs/400; b = Fs/80; % xac dinh khoang do tre 
    L = length(lag); 
    ViTriBatDauDuyet = 1;  
    % xac dinh vi tri bat dau duyet (o giua cua do thi tu tuong quan) 
    % xac dinh F0 ---------------------------------------------------- 
    YttqMax = -1; 
    for i = ViTriBatDauDuyet:L     
        if Yttq(i)>YttqMax && a<=lag(i) && lag(i)<=b      
            YttqMax = Yttq(i);     
            lagmax = i;     
            F0 = Fs/lag(lagmax); 
            % tinh tan so co ban F0     
        end
    end
    % ------------------------------------------------------------
end

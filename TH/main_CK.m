close all;

global K N;
% Parameters
K = 3;
N = 39;

%khoi tao phan 2
%Khung a:
file01_a = '..\TinHieuHuanLuyen\NguyenAmA\01_a.wav';
file02_a = '..\TinHieuHuanLuyen\NguyenAmA\02_a.wav';
file03_a = '..\TinHieuHuanLuyen\NguyenAmA\03_a.wav';
file04_a = '..\TinHieuHuanLuyen\NguyenAmA\04_a.wav';
file05_a = '..\TinHieuHuanLuyen\NguyenAmA\05_a.wav';
file06_a = '..\TinHieuHuanLuyen\NguyenAmA\06_a.wav';
file07_a = '..\TinHieuHuanLuyen\NguyenAmA\07_a.wav';
file08_a = '..\TinHieuHuanLuyen\NguyenAmA\08_a.wav';
file09_a = '..\TinHieuHuanLuyen\NguyenAmA\09_a.wav';
file10_a = '..\TinHieuHuanLuyen\NguyenAmA\10_a.wav';
file11_a = '..\TinHieuHuanLuyen\NguyenAmA\11_a.wav';
file12_a = '..\TinHieuHuanLuyen\NguyenAmA\12_a.wav';
file14_a = '..\TinHieuHuanLuyen\NguyenAmA\14_a.wav';
file15_a = '..\TinHieuHuanLuyen\NguyenAmA\15_a.wav';
file16_a = '..\TinHieuHuanLuyen\NguyenAmA\16_a.wav';
file17_a = '..\TinHieuHuanLuyen\NguyenAmA\17_a.wav';
file18_a = '..\TinHieuHuanLuyen\NguyenAmA\18_a.wav';
file19_a = '..\TinHieuHuanLuyen\NguyenAmA\19_a.wav';
file20_a = '..\TinHieuHuanLuyen\NguyenAmA\20_a.wav';
file21_a = '..\TinHieuHuanLuyen\NguyenAmA\21_a.wav';
file22_a = '..\TinHieuHuanLuyen\NguyenAmA\22_a.wav';

%Nguyen am e
file01_e = '..\TinHieuHuanLuyen\NguyenAmE\01_e.wav';
file02_e = '..\TinHieuHuanLuyen\NguyenAmE\02_e.wav';
file03_e = '..\TinHieuHuanLuyen\NguyenAmE\03_e.wav';
file04_e = '..\TinHieuHuanLuyen\NguyenAmE\04_e.wav';
file05_e = '..\TinHieuHuanLuyen\NguyenAmE\05_e.wav';
file06_e = '..\TinHieuHuanLuyen\NguyenAmE\06_e.wav';
file07_e = '..\TinHieuHuanLuyen\NguyenAmE\07_e.wav';
file08_e = '..\TinHieuHuanLuyen\NguyenAmE\08_e.wav';
file09_e = '..\TinHieuHuanLuyen\NguyenAmE\09_e.wav';
file10_e = '..\TinHieuHuanLuyen\NguyenAmE\10_e.wav';
file11_e = '..\TinHieuHuanLuyen\NguyenAmE\11_e.wav';
file12_e = '..\TinHieuHuanLuyen\NguyenAmE\12_e.wav';
file14_e = '..\TinHieuHuanLuyen\NguyenAmE\14_e.wav';
file15_e = '..\TinHieuHuanLuyen\NguyenAmE\15_e.wav';
file16_e = '..\TinHieuHuanLuyen\NguyenAmE\16_e.wav';
file17_e = '..\TinHieuHuanLuyen\NguyenAmE\17_e.wav';
file18_e = '..\TinHieuHuanLuyen\NguyenAmE\18_e.wav';
file19_e = '..\TinHieuHuanLuyen\NguyenAmE\19_e.wav';
file20_e = '..\TinHieuHuanLuyen\NguyenAmE\20_e.wav';
file21_e = '..\TinHieuHuanLuyen\NguyenAmE\21_e.wav';
file22_e = '..\TinHieuHuanLuyen\NguyenAmE\22_e.wav';


%Nguyen am i
file01_i = '..\TinHieuHuanLuyen\NguyenAmI\01_i.wav';
file02_i = '..\TinHieuHuanLuyen\NguyenAmI\02_i.wav';
file03_i = '..\TinHieuHuanLuyen\NguyenAmI\03_i.wav';
file04_i = '..\TinHieuHuanLuyen\NguyenAmI\04_i.wav';
file05_i = '..\TinHieuHuanLuyen\NguyenAmI\05_i.wav';
file06_i = '..\TinHieuHuanLuyen\NguyenAmI\06_i.wav';
file07_i = '..\TinHieuHuanLuyen\NguyenAmI\07_i.wav';
file08_i = '..\TinHieuHuanLuyen\NguyenAmI\08_i.wav';
file09_i = '..\TinHieuHuanLuyen\NguyenAmI\09_i.wav';
file10_i = '..\TinHieuHuanLuyen\NguyenAmI\10_i.wav';
file11_i = '..\TinHieuHuanLuyen\NguyenAmI\11_i.wav';
file12_i = '..\TinHieuHuanLuyen\NguyenAmI\12_i.wav';
file14_i = '..\TinHieuHuanLuyen\NguyenAmI\14_i.wav';
file15_i = '..\TinHieuHuanLuyen\NguyenAmI\15_i.wav';
file16_i = '..\TinHieuHuanLuyen\NguyenAmI\16_i.wav';
file17_i = '..\TinHieuHuanLuyen\NguyenAmI\17_i.wav';
file18_i = '..\TinHieuHuanLuyen\NguyenAmI\18_i.wav';
file19_i = '..\TinHieuHuanLuyen\NguyenAmI\19_i.wav';
file20_i = '..\TinHieuHuanLuyen\NguyenAmI\20_i.wav';
file21_i = '..\TinHieuHuanLuyen\NguyenAmI\21_i.wav';
file22_i = '..\TinHieuHuanLuyen\NguyenAmI\22_i.wav';


%Nguyen am O
file01_o = '..\TinHieuHuanLuyen\NguyenAmO\01_o.wav';
file02_o = '..\TinHieuHuanLuyen\NguyenAmO\02_o.wav';
file03_o = '..\TinHieuHuanLuyen\NguyenAmO\03_o.wav';
file04_o = '..\TinHieuHuanLuyen\NguyenAmO\04_o.wav';
file05_o = '..\TinHieuHuanLuyen\NguyenAmO\05_o.wav';
file06_o = '..\TinHieuHuanLuyen\NguyenAmO\06_o.wav';
file07_o = '..\TinHieuHuanLuyen\NguyenAmO\07_o.wav';
file08_o = '..\TinHieuHuanLuyen\NguyenAmO\08_o.wav';
file09_o = '..\TinHieuHuanLuyen\NguyenAmO\09_o.wav';
file10_o = '..\TinHieuHuanLuyen\NguyenAmO\10_o.wav';
file11_o = '..\TinHieuHuanLuyen\NguyenAmO\11_o.wav';
file12_o = '..\TinHieuHuanLuyen\NguyenAmO\12_o.wav';
file14_o = '..\TinHieuHuanLuyen\NguyenAmO\14_o.wav';
file15_o = '..\TinHieuHuanLuyen\NguyenAmO\15_o.wav';
file16_o = '..\TinHieuHuanLuyen\NguyenAmO\16_o.wav';
file17_o = '..\TinHieuHuanLuyen\NguyenAmO\17_o.wav';
file18_o = '..\TinHieuHuanLuyen\NguyenAmO\18_o.wav';
file19_o = '..\TinHieuHuanLuyen\NguyenAmO\19_o.wav';
file20_o = '..\TinHieuHuanLuyen\NguyenAmO\20_o.wav';
file21_o = '..\TinHieuHuanLuyen\NguyenAmO\21_o.wav';
file22_o = '..\TinHieuHuanLuyen\NguyenAmO\22_o.wav';


%Nguyen am U
file01_u = '..\TinHieuHuanLuyen\NguyenAmU\01_u.wav';
file02_u = '..\TinHieuHuanLuyen\NguyenAmU\02_u.wav';
file03_u = '..\TinHieuHuanLuyen\NguyenAmU\03_u.wav';
file04_u = '..\TinHieuHuanLuyen\NguyenAmU\04_u.wav';
file05_u = '..\TinHieuHuanLuyen\NguyenAmU\05_u.wav';
file06_u = '..\TinHieuHuanLuyen\NguyenAmU\06_u.wav';
file07_u = '..\TinHieuHuanLuyen\NguyenAmU\07_u.wav';
file08_u = '..\TinHieuHuanLuyen\NguyenAmU\08_u.wav';
file09_u = '..\TinHieuHuanLuyen\NguyenAmU\09_u.wav';
file10_u = '..\TinHieuHuanLuyen\NguyenAmU\10_u.wav';
file11_u = '..\TinHieuHuanLuyen\NguyenAmU\11_u.wav';
file12_u = '..\TinHieuHuanLuyen\NguyenAmU\12_u.wav';
file14_u = '..\TinHieuHuanLuyen\NguyenAmU\14_u.wav';
file15_u = '..\TinHieuHuanLuyen\NguyenAmU\15_u.wav';
file16_u = '..\TinHieuHuanLuyen\NguyenAmU\16_u.wav';
file17_u = '..\TinHieuHuanLuyen\NguyenAmU\17_u.wav';
file18_u = '..\TinHieuHuanLuyen\NguyenAmU\18_u.wav';
file19_u = '..\TinHieuHuanLuyen\NguyenAmU\19_u.wav';
file20_u = '..\TinHieuHuanLuyen\NguyenAmU\20_u.wav';
file21_u = '..\TinHieuHuanLuyen\NguyenAmU\21_u.wav';
file22_u = '..\TinHieuHuanLuyen\NguyenAmU\22_u.wav';



array_a = {file01_a, file02_a,file03_a, file04_a, file05_a,file06_a,file07_a,file08_a,file09_a,file10_a, file11_a,file12_a,file14_a,file15_a,file16_a,file17_a,file18_a, file19_a ,file20_a,file21_a,file22_a};
array_e = {file01_e, file02_e,file03_e, file04_e, file05_e,file06_e,file07_e,file08_e,file09_e,file10_e, file11_e,file12_e,file14_e,file15_e,file16_e,file17_e,file18_e,file19_e,file20_e,file21_e,file22_e};
array_i = {file01_i, file02_i,file03_i, file04_i, file05_i,file06_i,file07_i,file08_i,file09_i,file10_i, file11_i,file12_i,file14_i,file15_i,file16_i,file17_i,file18_i,file19_i,file20_i,file21_i,file22_i};
array_o = {file01_o, file02_o,file03_o, file04_o, file05_o,file06_o,file07_o,file08_o,file09_o,file10_o, file11_o,file12_o,file14_o,file15_o,file16_o,file17_o,file18_o,file19_o,file20_o,file21_o,file22_o};
array_u = {file01_u, file02_u,file03_u, file04_u, file05_u,file06_u,file07_u,file08_u,file09_u,file10_u, file11_u,file12_u,file14_u,file15_u,file16_u,file17_u,file18_u,file19_u,file20_u,file21_u,file22_u};
Cal(array_a,array_e,array_i,array_o,array_u);
function cm = Cal(array_a,array_e,array_i,array_o,array_u)
    correct = [];
    numOfVowel = 105;
    for j = 1:10
        % Create feature vector database
        db = CreateDB(array_a,array_e,array_i,array_o,array_u);
        % Predict test signals using created database
        [predict, trueVal] = PredictAll("..\NguyenAmKiemThu-16k", db);
        % Plot confusion matrix
        cm = confusionmat(trueVal, predict);
        confusionchart(cm, ["/a/", "/e/", "/i/", "/o/", "/u/"]);
        correct(j) = 0;
        for i = 1:105
            if trueVal(i) == predict(i)
                correct(j) = correct(j) + 1;
            end
        end
        correct(j) = correct(j) / numOfVowel*100;
        disp("Lan "+j+" co ti le chinh xac(%) la "+correct(j));
    end    
    giatritrungbinh = mean(correct);
    saiso = std(correct);
end 
function [threshHold] = Compute_Threshold(STE_PowFrame_Matrix, Weight)
    [histSTE, x_STE] = hist(STE_PowFrame_Matrix, round(length(STE_PowFrame_Matrix)/0.5)); % T?n su?t xu?t hi?n ( hist STE ) giá tr? STE m?i frame
    % t?i các v? trí x_STE.
    % vecto histSTE : l?u t?n su?t xu?t hi?n ( s? l?n xu?t hi?n ) giá tr? STE.
    % c?a m?i frame ( STE_PowFrame_Matrix) t?i v? trí x_STE ( vecto ).
    maximaHistSTE1 = 0;
    maximaHistSTE2 = 0;
    maximaIndex1 = 0; % V? trí c?c ??i c?c b? th? 1
    maximaIndex2 = 0; % V? trí c?c ??i c?c b? th? 2
    %Tìm c?c ??i c?c b? th? nh?t và th? hai n?m cùng 1 frame
    for i = 2 : length(histSTE) - 1  % Duy?t k?t qu? ?? th? t?n su?t ( histSTE)
        previous = i - 1;
        next = i + 1;
        while(histSTE(i) == histSTE(next)) % Xét v? trí histSTE th? i và histSTE li?n k?
            next = next + 1;
        end
        if(histSTE(i) > histSTE(previous) && histSTE(i) > histSTE(next)) % Ki?m tra giá tr? t?i histSTE th? i so v?i giá tr? t?i histSTE tr??c và sau
            if(maximaIndex1 == 0)
                maximaHistSTE1 = histSTE(i);
                maximaIndex1 = i;
            else
                maximaHistSTE2 = histSTE(i);
                maximaIndex2 = i;
                break;
            end
        end
        i = next;
    end
    maximaHistSTE1 = x_STE(maximaIndex1); % K?t qu? giá tr? c?c ??i c?c b? th? nh?t
    maximaHistSTE2 = x_STE(maximaIndex2); % K?t qu? giá tr? c?c ??i c?c b? th? hai
    % B2: Áp d?ng công th?c : T = (W * M1 + M2) / (W + 1)
    threshHold = (Weight * maximaHistSTE1 + maximaHistSTE2) / (Weight + 1);
end


function [sig,Fs1] = Speech_Silence(fileName)
    [x1,Fs1] = audioread(fileName);
    t = 0:1/Fs1:(length(x1)-1)/Fs1; %time of sample

    frame_len = round(0.02*Fs1); % làm tròn v? s? nguyên g?n nh?t
    N = length(x1);
    num_frame = floor(N/frame_len); %làm tròn v? s? nguyên nh? nh?t

    %frame by frame with 0.02
    for i=1:num_frame
        frame = x1(frame_len*(i-1)+1:frame_len*i);
         x1_2 = frame.*frame;
         STE(i)= sum(x1_2);
    end
    time_STE = 0:0.02:(N-frame_len)/Fs1;
    % standard signal

    x1 = x1./max(abs(x1));
    STE = STE./max(STE);

    T = Compute_Threshold(STE, 100);
    % T=0.108;
    %find silence
     ste_nor = STE>T;
    out = [];
     num = 0.3/0.02+1;
     %draw Amplitude to discriminate
     u=1;i_time=1;

     while (u <= length(ste_nor))
         d=0;
         inte = u;
          while (ste_nor(u)==0 && u<length(ste_nor))
             d = d+1;
             u=u+1;
          end     
           if (d>=num)  
               if (u<length(ste_nor))
                     u = u-1;
               end
               out(i_time) = inte;
               i_time = i_time+1;
               out(i_time) = u;
               i_time = i_time +1;
          end
          u=u+1;   
     end
     if (length(out)>=4)
     dau = floor((out(2)-1)*0.02*Fs1);
     cuoi = floor((out(3) -1)*0.02*Fs1);
     sig = x1(dau:cuoi);
     else
         sig = x1;
     end
end

function [sig] = findStable(sig0)
    sig0_length = length(sig0);
    n_divide = 3;
    sig_first = floor(sig0_length/n_divide);
    sig_last = floor(sig0_length/n_divide * 2);
    sig = sig0(sig_first:sig_last);
end

function [mfcc_frame]= mffc(sig,fs)
    %matric cua vector dac trung cua tat ca cac khung cua doan tin hieu nguyen
    %am duoc danh dau
    global N;
    mfcc_frame = melcepst(sig, fs, 'E', N-1, floor(3*log(fs)), 0.03*fs, 0.01*fs);
end

%ham tich vector dac trung cho 1 nguyen am cua 1 nguoi noi
function [featuredVector] = CalFeaturedVector(matrix)
    featuredVector = mean(matrix);
end

%ham tich vector dac trung cho 1 nguyen am cua nhieu nguoi noi
function [kq] = CalFeaturedVectorForPeople(array)
    mfcc_kq=[];      
    for i =1 : length(array)
        [output,fs] = Speech_Silence(array{i});
        [sig] = findStable(output);
        mfccOneFile = mffc(sig,fs);
        featuredVector = CalFeaturedVector(mfccOneFile);
        mfcc_kq = [mfcc_kq;featuredVector];
    end
    kq = mean(mfcc_kq);
end


function [arrayCoeff] = matrixMfccTotalFrame(array)
    arrayCoeff = [];
     for i =1 : length(array)
        [output,fs] = Speech_Silence(array{i});
        sig = findStable(output);
        mfccOneFile = mffc(sig,fs);
        arrayCoeff = [arrayCoeff;mfccOneFile];
    end
end

function [mfcc_a,mfcc_e,mfcc_i,mfcc_o,mfcc_u] = matrixMfccTotalFrameInTotalFile(array_a,array_e,array_i,array_o,array_u)
    mfcc_a=matrixMfccTotalFrame(array_a);      
    mfcc_e=matrixMfccTotalFrame(array_e);      
    mfcc_i=matrixMfccTotalFrame(array_i);      
    mfcc_o=matrixMfccTotalFrame(array_o);      
    mfcc_u=matrixMfccTotalFrame(array_u);      
end

% Phat
% Create feature vector database (using kmeans in Statistics and Machine Learning Toolbox)
function db = CreateDB(array_a,array_e,array_i,array_o,array_u)
    global K;
    [mfcc_a,mfcc_e,mfcc_i,mfcc_o,mfcc_u] = matrixMfccTotalFrameInTotalFile(array_a,array_e,array_i,array_o,array_u);
    [idx_a,db(:,:,1),sumd_a,D_a] = kmeans(mfcc_a, K);
    [idx_e,db(:,:,2),sumd_e,D_e] = kmeans(mfcc_e, K);
    [idx_i,db(:,:,3),sumd_i,D_i] = kmeans(mfcc_i, K);
    [idx_o,db(:,:,4),sumd_o,D_o] = kmeans(mfcc_o, K);
    [idx_u,db(:,:,5),sumd_u,D_u] = kmeans(mfcc_u, K);
end

% Predict vowel (1 signal) based on given database
function [vowel,euc_dis] = Predict1(filename, db)
    global K;
    % Find mfcc vectors and their mean on 1 signal
    [output,fs] = Speech_Silence(filename);
    sig = findStable(output);
    current_mffc = mffc(sig,fs);
    current_fvector = CalFeaturedVector(current_mffc);
    
    euc_dis = zeros(K,5);
    for i=1:K
        for j=1:5
            % Calculate Euclidian distances
            euc_dis(i,j) = sqrt(sum((current_fvector - db(i,:,j)) .^ 2));
        end
    end
    
    [m,idx] = min(euc_dis,[],2); % min of rows
    [m_fin, idx_fin] = min(m);   % min of columns (overall min)
    vowel = idx(idx_fin);
end

% Predict all test signals 
function [predict,trueVal] = PredictAll(testfolder, db)
    % Read all files in testfolder
    % and return their true and predicted vowels (in numbers 1-5 for a,e,i,o,u respectively)
    dir_content = dir(testfolder);
    count = 0;
    trueVal = [];
    predict = [];
    for i=3:length(dir_content)
        subdir_name = getfield(dir_content(i), "name");
        subdir_content = dir(testfolder + "\" + subdir_name);
        for j=3:length(subdir_content)
            filename = getfield(subdir_content(j), "name");
            if filename == "a.wav"
                count = count + 1;
                trueVal(count) = 1;
                predict(count) = Predict1(testfolder + "\" + subdir_name + "\" + filename, db);
                continue
            end
            if filename == "e.wav"
                count = count + 1;
                trueVal(count) = 2;
                predict(count) = Predict1(testfolder + "\" + subdir_name + "\" + filename, db);
                continue
            end
            if filename == "i.wav"
                count = count + 1;
                trueVal(count) = 3;
                predict(count) = Predict1(testfolder + "\" + subdir_name + "\" + filename, db);
                continue
            end
            if filename == "o.wav"
                count = count + 1;
                trueVal(count) = 4;
                predict(count) = Predict1(testfolder + "\" + subdir_name + "\" + filename, db);
                continue
            end
            if filename == "u.wav"
                count = count + 1;
                trueVal(count) = 5;
                predict(count) = Predict1(testfolder + "\" + subdir_name + "\" + filename, db);
            end
        end
    end
end
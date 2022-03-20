% Bien chuan trong file .lab
t44=[0 0.93 1.42 2.59 3 4.71 5.11 6.26 6.66 8.04 8.39 9.27];
t30 =[0 0.59 0.97 1.76 2.11 3.44 3.77 4.7 5.13 5.96 6.28 6.78];
t42 = [0 0.46 0.99 1.56 2.13 2.51 2.93 3.79 4.38 4.77 5.22 5.79];
t45 = [0 0.88 1.34 2.35 2.82 3.76 4.13 5.04 5.5 6.41 6.79 7.42];

t06=[1.52 1.92 3.91 4.35 6.18 6.6 8.67 9.14 10.94 11.33];
t01=[0.45 0.81 1.53 1.85 2.69 2.86 3.78 4.15 4.84 5.14];
t02=[0.83 1.37 2.09 2.6 3.57 4 4.76 5.33 6.18 6.68];
t03=[1.03 1.42 2.46 2.8 4.21 4.52 6.81 7.14 8.22 8.5];

file_01 = '01MDA.wav';
file_02 = '02FVA.wav';
file_03 = '03MAB.wav';
file_06 = '06FTB.wav';
file_30 = '30FTN.wav';
file_42 = '42FQT.wav';
file_44 = '44MTT.wav';
file_45 = '45MDV.wav';

% Thuc hien 4 file tin hieu 
[mark,t_sample2,STE,mark_silence] = Discriminatory(file_01,t01);
figure('Name',file_01);
findFrequency(file_01,mark,t_sample2,STE,mark_silence,t01);
[mark,t_sample2,STE,mark_silence] = Discriminatory(file_02,t02);
figure('Name',file_02);
findFrequency(file_02,mark,t_sample2,STE,mark_silence,t02);

[mark,t_sample2,STE,mark_silence]= Discriminatory(file_03,t03);
figure('Name',file_03);
findFrequency(file_03,mark,t_sample2,STE,mark_silence,t03);


[mark,t_sample2,STE,mark_silence]=Discriminatory(file_06,t06);
figure('Name',file_06);
findFrequency(file_06,mark,t_sample2,STE,mark_silence,t06);


[mark,t_sample2,STE,mark_silence] = Discriminatory(file_44,t44);
figure('Name',file_44);
findFrequency(file_44,mark,t_sample2,STE,mark_silence,t44);

[mark,t_sample2,STE,mark_silence] = Discriminatory(file_45,t45);
figure('Name',file_45);
findFrequency(file_45,mark,t_sample2,STE,mark_silence,t45);

[mark,t_sample2,STE,mark_silence]= Discriminatory(file_30,t30);
figure('Name',file_30);
findFrequency(file_30,mark,t_sample2,STE,mark_silence,t30);


[mark,t_sample2,STE,mark_silence]=Discriminatory(file_42,t42);
figure('Name',file_42);
findFrequency(file_42,mark,t_sample2,STE,mark_silence,t42);






%ham phan biet nguyen am va khoang lang
function [mark,t_sample2,STE,mark_silence] = Discriminatory(sig,t_lab)
% doc file tin hieu, khoi tao cac gia tri ban dau
[x,fs] = audioread(sig);
tinite = 0.03;
half_tinite=tinite/2;
time_shift= 0.01;
m_sample = ceil(tinite*fs);
time_sig = length(x)/fs; 
amount_fragment = round(length(x)/(m_sample/3))-2;
time_frame_first =0;
time_frame_last = 0.03;
T=0.006;% T la nguong

%tao ma tran khung tin hieu
for i = 1:amount_fragment-1
    sample_first = time_frame_first*fs+1;
    sample_last = time_frame_last*fs;
    n_test = x(sample_first:sample_last);
    for j=1:m_sample
    n_frame(i,j)=n_test(j);
    end
    time_frame_first = time_frame_last-0.02;
    time_frame_last = time_frame_first+0.03; 
end


% xu li khung tin hieu cuoi cung
sample_first = time_frame_first*fs+1;
n_test = x(sample_first:length(x));
for j=1:length(n_test)
    n_frame(amount_fragment,j)=n_test(j);
end 
for j=length(n_test)+1:m_sample
    n_frame(amount_fragment,j)=0;
end 


%STE
[rows,cols] = size(n_frame);
for i =1:rows
      STE(i) = sum(n_frame(i,:).^2)
      
end
% chuan hoa STE ve [0-1]
STE = STE/max(STE);
P= STE>T;
kk=0;
i=1;
% danh dau khoang lang
while (i<=length(P))
    if (P(i)==0)
        index_dau = i;
        for j=i:length(P)
            if (P(j)==1)
                break;
            end
        end
        index_cuoi =j-1;
        time_dau = (index_dau-1)*time_shift;
        time_cuoi = 2*half_tinite + (index_cuoi-1)*time_shift;
        time_silence = time_cuoi-time_dau;
        if (time_silence>0.3)
            kk=kk+1;
            mark_silence(kk)=time_dau;
            kk=kk+1;
            mark_silence(kk)=time_cuoi;
        end
        i=j;
        
    end
    i=i+1;
end
% danh dau nguyen am
for i=2:length(mark_silence)-1
    mark(i-1)=mark_silence(i)+time_shift;
end


t_sample1 = 0.015:0.01:time_sig;
if (t_sample1(length(t_sample1))+0.015>time_sig)
    t_sample2 = t_sample1(1:length(t_sample1)-1);
end
end

%Tim tan so co ban 
function findFrequency(file,mark,t_sample2,STE,mark_silence,t_lab)
% doc file tin hieu 
[sig,fs] = audioread(file); 
N = 32768; 
frame_len = ceil(0.03*fs); 
% ham cua so hamming
h = hamming(frame_len);
total_frame =round(length(sig)/(frame_len/3))-2;
time_frame_first =0;
time_frame_last = 0.03;
s=0;k=0;
for i = 1 : total_frame-1
    %kiem tra khung tin hieu co phai la nguyen am khong, ktra dung thi tim
    %f0, sai thi bo qua
   is_speech =checkSpeech(mark,time_frame_first,time_frame_last);
   if (is_speech==1)
    time_output(i)=(time_frame_first+time_frame_last)/2;
    sample_first = time_frame_first*fs+1;
    sample_last = time_frame_last*fs;
    % lay mau trong 1 khung tin hieu
    sample_frame = sig(sample_first:sample_last);
    frame = h.*sample_frame;
    %ham fft
    P2 = abs(fft(frame,N));
    P1 = P2(1:length(P2)/2+1); 
    % downsample
    data2 = downsample_signal(P1,2);
    data3 = downsample_signal(P1,3);
    data4 = downsample_signal(P1,4);
    data5 = downsample_signal(P1,5);
    data6 = downsample_signal(P1,6);
    data7 = downsample_signal(P1,7);
    data8 = downsample_signal(P1,8);
    data9 = downsample_signal(P1,9);
    
    %kiem tra dieu kien de xac dinh so lan nen mau
  if (checkAmplitude(P1)==1)
    for j = 1:length(data9)
     hps(j) = P1(j)*data2(j)*data3(j)*data4(j)*data5(j)*data6(j)*data7(j)*data8(j)*data9(j);
    end
      
  else
      for j = 1:length(data3)
      hps(j) = P1(j)*data2(j)*data3(j);
      end
  
  end
   % tim gia tri tan so co bien do lon nhat  
     [m,n]=findpeaks(hps, 'SORTSTR', 'descend');
     F0 =  ( (n(1) / 32768) * fs );
  
   % kiem tra F0 thoa man 70<F0<400
     if (F0>70 && F0<400)
     output(i) = F0;
     k=k+1;
     output1(k)=F0;
     else
     output(i)=0;
     end
   
   else
           time_output(i)=(time_frame_first+time_frame_last)/2;
           output(i)=0;
   end
   time_frame_first = time_frame_last-0.02;
   time_frame_last = time_frame_first+0.03;
end
%mean va std
[mean,std]=calculation_F0(output);
%plot
t_sample = 0:1/fs:(length(sig)-1)/fs;
%do thi phan biet nguyen am va khoang lang bao gom tin hieu dau vao, STE,
%bien chuan va bien tim duoc nho thuat toan
subplot(4,1,1);
plot(t_sample,sig); title('Speech/Silence'); hold on;
plot(t_sample2,STE,'r','LineWidth',2);hold on;
for i =1:length(t_lab)
    plot ([t_lab(i) t_lab(i)],[-1 1],'r');
    hold on;
end
for i =1:length(mark_silence)
    plot ([mark_silence(i) mark_silence(i)],[-1 1],'b');
    hold on;
end

for i =1:length(mark)
    mark1(i)=mark(i)*100;
end
xlabel('time');ylabel('Magnitude');
legend('Algorithm' , 'Standard');


%do thi HPS cho 1 khung tin hieu nguyen am cuoi cung duoc xet
subplot(4,1,2)
 k=1:N;
 w=k*fs/N; % frequency axis
plot(w(1:1000),hps(1:1000));xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('HPS');

% do thi F0
subplot(4,1,3)
for i =1:length(mark1)
    plot ([mark1(i) mark1(i)],[0 400],'b');
    hold on;
end
plot(output,'.');
title(['F0 mean = ',num2str(mean),'F0 std =',num2str(std)] );
xlabel('Sample');
ylabel('Magnitude');


subplot(4,1,4)
fft_F0 = MedSmoothing(output,7);
[mean1,std1]=calculation_F0(fft_F0);

plot(fft_F0,'.');
title(['F0 mean = ',num2str(mean1),'F0 std =',num2str(std1)] );

end
%Ham kiem tra khung tin hieu co phai nguyen am
function [conclude] = checkSpeech(arr_speech,frame_time_first,frame_time_last)
    conclude=0;
    for i = 1:length(arr_speech)
       if (mod(i,2)==1 && (frame_time_first>=arr_speech(i) && frame_time_last<=arr_speech(i+1)))
           conclude = 1;
           break;
       end
    end

end

 %Ham kiem tra cac tan so co bien do cao trong 1 khoang
 function [conclude] = checkAmplitude(fft_signal)
    [m1,n1]=findpeaks(fft_signal, 'SORTSTR', 'descend');
    n2=n1(1:5);
    conclude=0;
    for i = 1:length(n2)
        if (n2(i)>740)
            conclude=1;
            break;
        end
    end
 end
 
 % ham downsample
 function [array] = downsample_signal(input,n)
    i=1;
    k=0;
    while (i<=length(input))
        k=k+1;
        array(k)=input(i);
        i=i+n;
    end
 end
 
 %ham tinh mean va std c?a F0
 function [mean,std] = calculation_F0(input)
    k=0;sum=0;std=0;
    for i=1:length(input)
        if (input(i)>0)
            sum=sum+input(i);
            k=k+1;
        end
    end
    mean = sum/k;
    for i=1:length(input)
          if (input(i)>0)
        std=std+(mean-input(i))*(mean-input(i));
          end
    end
    std=std/k;
    std = sqrt(std);
 end
 
 function y = MedSmoothing(x,N) 
 y = x;
 temp = 1:N; 
 for k = 1:length(x)
 for i = 1:N 
 if(k < ceil(N/2))
 if(i < ceil(N/2)-k+1) 
 temp(i) = 0;
 else
 temp(i) = x(k+i-ceil(N/2));
 end
 else
 if(k > length(x) - ceil(N/2) + 1)
 if(i < length(x) + ceil(N/2) - k + 1)
 temp(i) = x(k-ceil(N/2) + i);
 else
 temp(i) = 0; 
 end
 else
 temp(i) = x(k-ceil(N/2)+i); 
 end
 end
 end
 y(k) = median(temp);
 end
 end


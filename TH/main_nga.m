close all;
close all;

%audioread('30FTN.wav');
fileName1= '30FTN.wav';
FTN = [0.00 0.59 0.97 1.76 2.11 3.44 3.77 4.70 5.13 5.96 6.28 6.78];

%audioread('42FQT.wav');
fileName2= '42FQT.wav';
FQT = [0.00 0.46 0.99 1.56 2.13 2.51 2.93 3.79 4.38 4.77 5.22 5.79];

%audioread('44MTT.wav');
fileName3= '44MTT.wav';
MTT = [0.00 0.93 1.42 2.59 3.00 4.71 5.11 6.26 6.66 8.04 8.39 9.27];

%audioread('45MDV.wav');
fileName4= '45MDV.wav';
MDV = [0.00 0.88 1.34 2.35 2.82 3.76 4.13 5.04 5.50 6.41 6.79 7.42];

%audioread('01MDA.wav')
fileName5 = '01MDA.wav';
MDA = [0.00 0.45 0.81 1.53 1.85	2.69 2.86 3.78 4.15 4.84 5.14 5.58];

%audioread('02FVA.wav')
fileName6 = '02FVA.wav';
FVA = [0.00	0.83 1.37 2.09 2.60	3.57 4.00 4.76 5.33 6.18 6.68 7.18];

%audioread('03MAB.wav')
fileName7 = '03MAB.wav';
MAB = [0.00	1.03 1.42 2.46 2.80	4.21 4.52 6.81 7.14	8.22 8.50 9.37];

%audioread('06FTB.wav')
fileName8 = '06FTB.wav';
FTB = [0.00	1.52 1.92 3.91 4.35 6.18 6.60 8.67 9.14	10.94 11.33	12.75];

%draw
[x1,Fs1] = audioread(fileName1);
[t,time_STE,STE,time_mark]=Speech_Silence(fileName1,FTN);
figure(1);
caub(fileName1,FTN,t,time_STE,STE,time_mark);


[x2,Fs2] = audioread(fileName2);
[t,time_STE,STE,time_mark]=Speech_Silence(fileName2,FQT);
figure(2);
 caub(fileName2,FQT,t,time_STE,STE,time_mark);



[x3,Fs3] = audioread(fileName3);
[t,time_STE,STE,time_mark]=Speech_Silence(fileName3,MTT);
figure(3);
caub(fileName3,MTT,t,time_STE,STE,time_mark);


[x4,Fs4] = audioread(fileName4);
[t,time_STE,STE,time_mark]=Speech_Silence(fileName4,MDV);
figure(4);
caub(fileName4,MDV,t,time_STE,STE,time_mark);

[x5,Fs5] = audioread(fileName5);
[t,time_STE,STE,time_mark]=Speech_Silence(fileName5,MDA);
figure(5);
caub(fileName5,MDA,t,time_STE,STE,time_mark);

[x6,Fs6] = audioread(fileName6);
[t,time_STE,STE,time_mark]=Speech_Silence(fileName6,FVA);
figure(6);
caub(fileName6,FVA,t,time_STE,STE,time_mark);

[x7,Fs7] = audioread(fileName7);
[t,time_STE,STE,time_mark]=Speech_Silence(fileName7,MAB);
figure(7);
caub(fileName7,MAB,t,time_STE,STE,time_mark);

[x8,Fs8] = audioread(fileName8);
[t,time_STE,STE,time_mark]=Speech_Silence(fileName8,FTB);
figure(8);
caub(fileName8,FTB,t,time_STE,STE,time_mark);


function [t,time_STE,STE,time_mark] = Speech_Silence(fileName,standard)
[x1,Fs1] = audioread(fileName);
t = 0:1/Fs1:(length(x1)-1)/Fs1; %time of sample
frame_len = round(0.02*Fs1); % l�m tr�n v? s? nguy�n g?n nh?t
half = round(frame_len/2); %For overlapping frame 
N = length(x1);

%frame by frame with 0.02
for k = 1 : length(x1)/half -1 
    range = (k-1)*half + 1:(frame_len + (k-1)*half);
    frame = x1(range);
    x1_2 = frame.*frame;
    STE(k)= sum(x1_2);
end
time_STE = 0:0.01:(N-frame_len)/Fs1;
% standard signal

x1 = x1./max(abs(x1));
STE = STE./max(STE);

 T=0.0308 ;
%find silence
 ste_nor = STE>T;

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
       if (d>num)  
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
 for i =1:length(out)
     time_mark(i) = (out(i)-1)*0.01;
 end
 
  
end


function  caub(filename,standard,t,time_STE,STE,time_mark)
[sig,fs] = audioread(filename); %doc tin hieu
N = 2^15; %chon N co dinh
frame_len = round(0.03*fs); %do dai cua frame (30ms)
half = round(frame_len/2);
h = hamming(frame_len); %ham cua so hamming
i=1; %vi tri cua cac phan tu trong yyy_F0
PUnvoice = [];
for k = 1 : length(sig)/half -1 %Vong lap cac frame
    yyy_F0(k)=0;
 range = (k-1)*half + 1:(frame_len + (k-1)*half); %vi tri cac vong lap trong cua so
 frame = h.*sig(range); %gia tri cua moi phan tu trong window

 %dung ham FFT de phan tich pho cua frame
 P2 = abs(fft(frame,N)); 
 P1 = P2(1:length(P2)/2+1); 
  
 freq=linspace(0,fs/2,length(P1));
[y_value,y_peak] = findpeaks(P1,freq,'MinPeakHeight',2);
 if (length(y_peak)>=3)
 z_a = y_peak(3) - y_peak(2); 
 z_b = y_peak(2) - y_peak(1);
 if z_a<400 && z_a>70 %70Hz<F0<400Hz
 if z_b<400 && z_b>70
 yyy_F0(k) = (z_a+z_b)/2;%ket qua cua f0
 i=i+1;
 
 end
 end
 end
end

dau = standard(3);
cuoi = dau+0.03;
s_dau = dau*44100+1;
s_cuoi = cuoi*44100;
frame = h.*sig(s_dau:s_cuoi); %gia tri cua moi phan tu trong window
 %dung ham FFT de phan tich pho cua frame
P2 = abs(fft(frame,N)); 
P1 = P2(1:length(P2)/2+1); 

     PUnvoice = P2(1:length(P2)/2+1);    
     fregUnvoice=linspace(0,fs/2,length(PUnvoice));
 
%
dau = standard(2);
cuoi = dau+0.03;
s_dau = dau*44100+1;
s_cuoi = cuoi*44100;
frame = h.*sig(s_dau:s_cuoi); %gia tri cua moi phan tu trong window
 %dung ham FFT de phan tich pho cua frame
P2 = abs(fft(frame,N)); 
P1 = P2(1:length(P2)/2+1); 

   Pvoice = P2(1:length(P2)/2+1);    
     fregvoice=linspace(0,fs/2,length(Pvoice));


%
%

subplot(4,2,[1,2]);
plot(time_STE,STE,'Color','red');
title('STE');
subplot(4,2,[3,4]);
plot(t,sig);
hold on;
plot(time_STE,STE,'Color','red');
hold on;

  for i =1:length(time_mark)
      plot ([time_mark(i) time_mark(i)],[-1 1],'b','linewidth',1.4);
      hold on;  
  end
   hold on;
    for i=1:length(standard)
       plot([standard(i) standard(i)],[-1 1],'r','linewidth',1.4);
       hold on;
    end
   legend('Algorithm' , 'Standard');
   
subplot(4,2,5);
plot(fregvoice(1:length(fregvoice)), Pvoice(1:length(Pvoice))); %V� t?n s? Fs = 441000Hz , gi?i h?n < 2000Hz ,
%m� t�nh peak tr�n n?a ph? .
% 44100 / 2 = 22050 . Sau ?� chia ti?p cho 10 ?? c� gi� tr? < 2000Hz
findpeaks(P1(1:length(P1)), freq(1:length(freq)), 'MinPeakHeight', 2);
title('Frame Voice');

subplot(4,2,6);
plot(fregUnvoice(1:length(fregUnvoice)/10), PUnvoice(1:length(PUnvoice)/10));
title('Frame Unvoice');

subplot(4,2,[7,8]);
plot(yyy_F0,'.');
 title('F0');
end
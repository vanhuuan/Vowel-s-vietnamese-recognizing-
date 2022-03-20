% chuong trinh lam tron 1 tin hieu de khu nhieu (denoise)
% dung bo loc trung binh truot 3 diem (3-points moving-averaging filter)
% co PTSP: y[n] = 1/3*(x[n-1]+x[n]+x[n+1]) (he ko nhan qua) 
%-> h[n]=1/3*(delta[n-1]+delta[n]+delta[n+1]) = [1/3, 1/3, 1/3] (n=-1,0,1)
%
% problem statement: co tin hieu bi lan voi nhieu x[n], di khoi phuc lai
% tin hieu goc s[n] (an trong x[n])
% ky vong KQ: y[n] cang giong s[n] cang tot

% sinh tin hieu bi lan voi nhieu cong (additive noise)
clear all;
clf;                            % clear figures
A = 0.5;                    % cong suat nhieu ti le voi sqr(A)  
L = 51;                         % do dai tin hieu
n = 0:L-1;                      % bien thoi gian roi rac
d = A*randn(1,L);             % sinh tin hieu Gaussian white noise d[n] (A=std cua tin hieu nhieu)
s = 2*n.*(0.9.^n);              % sinh tin hieu goc s[n] = 2n(0.9)^n
x = s + d;                      % tin hieu co nhieu x[n]=s[n]+d[n]

figure(1)                    
hold on
subplot(3,1,1)                 
plot(n,d,'r-',n,s,'k--',n,x,'b-.'); % ve do thi d[n],s[n],x[n]
xlabel('Chi so thoi gian n');
ylabel('Bien do');
legend('d[n]','s[n]','x[n]');
title('Noise d[n] vs. original s[n] vs. noisy signals x[n]');

% cach 1: cai dat he thong bang cach dich thoi gian va tinh
% TBC cua 3 tin hieu
x1 = [x(2:L), 0];     % x1[n] = x[n+1]
x2 = [x] ;                 % x2[n] = x[n]
x3 = [0, x(1:L-1)];    % x3[n] = x[n-1]

subplot(3,1,2)
holdd on
% ve do thi x[n-1],x[n],x[n+1]
plot(n,x1,'r-.',n,x2,'b-.',n,x3,'k-.');
xlabel('Chi so thoi gian n');
ylabel('Bien do');
legend('x[n+1]','x[n]','x[n-1]');
title('time-shifted signals of x[n]');

y1 = 1/3*(x1+x2+x3);     % lenh cai dat he thong
subplot(3,1,3)
hold on
% ve do thi y1[n] vs. s[n]
plot(n,y1(1:L),'r-',n,s(1:L),'b-');
xlabel('Chi so thoi gian n');
ylabel('Bien do');
legend('y1[n]','s[n]');
title('3-points smoothed y1[n] vs. original signal s[n]');

% cach 2: cai dat he thong bang cach dung ham tinh tong chap conv()
% he lay TB truot KHONG nhan qua = ghep noi tiep he som 1 don vi co h1[n]=delta(n+1) 
% voi he lay TB truot nhan qua co PTSP yNQ[n] = 1/3*(x[n]+x[n+1]+x[n+2])
hNQ = 1/3 * ones(1,3);    % hNQ[n] = [1/3, 1/3, 1/3] (Matlab hieu la n=0,1,2)
y2 = conv(x1, hNQ);         % y2[n] = x[n] * h[n] = x[n] * (h1[n] * hNQ[n]) = (x[n]*h1[n]) * hNQ[n] = x[n+1] * hNQ[n] 
% ve do thi y2[n] vs. s[n]
figure(2)
hold on
plot(n,y2(1:L),'r-',n,s(1:L),'b-');
xlabel('Chi so thoi gian n');
ylabel('Bien do');
legend('y2[n]','s[n]');
title('3-points smoothed y2[n] vs. original signal s[n]');

% ve do thi xep chong y1[n] va y2[n] de test ket qua
figure(3)
hold on
plot(n,y1,'r-',n,y2(1:L),'b-');
xlabel('Chi so thoi gian n');
ylabel('Bien do');
legend('y1[n]','y2[n]');
title('cach1  vs. cach2');

% cach 3: dung vong lap (hieu n la mot chi so mau nao do)
y3 = ones(1,L);
for i = 1:L
    if i == 1
        y_new = x(i)+x(i+1);
        y3(i) = 1/3*y_new;
        continue;
    end
    if i == L
        y_new = x(i)+x(i-1);
        y3(i) = 1/3*y_new;
        continue;
    end
    y_new = x(i-1)+x(i)+x(i+1);
    y3(i) = 1/3*y_new;
end
disp(x);
figure(4);
hold on
subplot(1,1,1)
% ve do thi y3[n] vs. s[n]
plot(n,y3(1:L),'r-',n,s(1:L),'b-');
xlabel('Chi so thoi gian n');
ylabel('Bien do');
legend('y3[n]','s[n]');
title('3-points smoothed y3[n] vs. original signal s[n]');

% ve do thi xep chong y1[n] , y2[n] va y3[n] de test ket qua
figure(5)
hold on
plot(n,y1,'r-',n,y2(1:L),'b-',n,y3(1:L),'c-');
xlabel('Chi so thoi gian n');
ylabel('Bien do');
legend('y1[n]','y2[n]','y3[n]');
title('cach1  vs cach2 vs cach 3');




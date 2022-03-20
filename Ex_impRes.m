% Cho he thong TTBB co PTSP: y(n) + 0.9y(n - 1) - 0.9y(n - 2) = 0.3x(n) + 0.2x(n - 1) - 0.3x(n - 2) 
% a. Tim DUXDV h[n] cua he
clf;
N = 1000;
num = [0.3 0.2 -0.3];     % numerator vector = {bj} = cac hs lien quan den x[n]
den = [1 0.9 -0.9];         % denominator vector = {ai} = cac hs lien quan den y[n] (a0=1)
h = impz(num,den,N);    % N: so luong mau cua h[n] (h[n] co the dai vo han)
%m = 0:N;
%h =  0.5.^m;    %h[n] = 0.5^m.u[m]: he on dinh
figure(1)
stem(h, 'fill');
title('h[n]');

% b. Tim y[n] cua he khi dc kich thich boi 1 x[n] nao do
n = 1:10^6;
x = randn(1,length(n));  % sinh tin hieu ngau nhien theo phan bo chuan co do dai 150 mau
%x = 10.^n;                    % sinh tin hieu co bien do kha lon co do dai 150 mau
y = conv(x,h);               % tinh convolution
figure(2)
hold on
subplot(311),stem(x(end-100:end),'fill');title('x[n]');
subplot(312),stem(h(end-100:end),'fill');title('h[n]');
subplot(313),stem(y(end-100:end),'fill');title('y[n]');
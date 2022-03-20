%   Bai 4

% Sine Parameters 
% x = A*cos(2*pi*F0*t+phase)
A = 1
F0 = 1000
phase = 0

% Sampling
stop = 2
%% x1[n] : Can recover
Fs1 = 3*F0
Ts1 = 1/Fs1
n1 = 0:Ts1:stop-Ts1
x1 = A*cos(2*pi*F0*n1+phase)
plot(n1(1:100),x1(1:100))
sound(x1,Fs1)
%% x2[n] : Cannot recover
Fs2 = 1.5*F0
Ts2 = 1/Fs2
n2 = 0:Ts2:stop-Ts2
x2 = A*cos(2*pi*F0*n2+phase)
plot(n2(1:100),x2(1:100))
sound(x2,Fs2)



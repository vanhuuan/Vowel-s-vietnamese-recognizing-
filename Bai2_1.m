% Bai 1
%% x[n]
n = -4:6
x = [0,0,0,-2,4,-2,2,-1,0,0,0]
stem(n,x,"filled")

%% 1a)  3x[-n-1]
n1 = -1-n
x1 = 3*x
stem(n1,x1,"filled")

%% 1b)  x[2n]-1
n2 = n/2
x2 = x-1
stem(n2,x2,"filled")

%% 1c) -x[n]+2
x3 = -x+2
stem(n,x3,"filled")
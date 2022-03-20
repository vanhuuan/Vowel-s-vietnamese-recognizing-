t = -10:10;  %%Can change the interval time by replacing 1 with 0.1
x = t>=0  %% For x[n]
xa = t>=2 % x[n-2]
xb = t>=7 % x[n-7]
xc = t<=0 % x[-n]
x1=2*xa
%stem(t,xc)    %%scatter can be used instead of plot
xlabel('time');
ylabel('amplitude');
title('x(n)=u(n)-u(n-3)');
 n = -5:1:20;
 xn = ((4/5).^n).*heaviside(n);
 stem(n,xn)
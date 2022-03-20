% Bai 2

%% y[n] = 2u[n-2] - 2u[n-7] - 2u[-n] + 2u[4-n]       (-15<=n<=15)
n = [-10:10]
y = 2*uss(2,-10,10) - 2*uss(7,-10,10) - 2*uss1(0,-10,10) + 2*uss1(-4,-10,10)
stem(n,y,"filled")

%% 2a)  2-3y[n]
x = 2-3*y
stem(n,x,"filled")

%% 2b)  3y[n-2]
n1 = n+2
x1 = 3*y
stem(n1,x1,"filled")

%% 2c)  2-2y[-2+n]
n2 = n+2
x2 = 2-2*y
stem(n2,x2,"filled")

%% Unit step sequence
% u[n-n0] (n1<=n,n0<=n2)
function [u,n] = uss(n0,n1,n2)
n = [n1:n2];
u = [n>=n0];
end

% u[-n-n0] (n1<=n,n0<=n2)
function [u,n] = uss1(n0,n1,n2)
n = [n1:n2];
u = [n<=-n0];
end
n = -4:6;
x = [0,0,0,-2,4,-2,2,-1,0,0,0];
tiledlayout(2,1)
nexttile
stem(n,x,"fill")
title('x')
nexttile
na = -n-1;
xa = 3*x;
stem(na,xa,"fill")
title('3x[-n-1]')
nb = n*2;
xb = x-1;
stem(nb,xb,"fill")
title('x[2n]-1')
nc = n;
xc = -x+2;
stem(nc,xc,"fill")
title('-x[n]+2')

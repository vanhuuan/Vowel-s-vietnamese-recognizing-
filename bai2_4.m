 f0 = 1000;
 fs1 = f0*3;
 fs2 = f0*1.5;
 fs = 512;  
 dt1 = 1/fs1; 
 dt2 = 1/fs2;
 StopTime = 2; 
 t1 = (0:dt1:StopTime);
 t2 = (0:dt2:StopTime);
 x1 = cos(2*pi*f0*t1);
 x2 = cos(2*pi*f0*t2);
 tiledlayout(2,1)
 nexttile
 plot(t1(1:100),x1(1:100),".-")
 title('x1[n]')
 nexttile
 plot(t2(1:100),x2(1:100),".-")
 title('x2[n]')
 disp("x1")
 sound(x1,fs1);
 pause(length(x1)/fs1);
 disp("x2")
 sound(x2,fs2);
 

[x, Fs] = audioread('lab_male.wav');
info = audioinfo('lab_male.wav');
t = (0:0.02:info.Duration); %Duration la lay thoi gian cua tin hieu
t = [t mod(info.Duration,0.02)];
t_x = 0:(1/Fs):info.Duration;
K = length(t) - 1; % so khung
arr = zeros(1,(length(x)));
%n = 0: N-1;
%x = sin(2*pi*n/Fs);
N = 0.02*Fs; %So mau trong 1 khung
for k =1 :K
    RXX = autom1(x, k, Fs, K);
    if(length(RXX) > 1) %kiểm tra khung cuoi cung co du hay ko
        r0 = RXX(1);
        nguong = 0.3*r0; %nguong thuong bang 30% cuc dai toan cuc dau tien
        peak = [];
        %location = [];
        for t = 2 : length(RXX)-1
            if(RXX(t)>RXX(t-1) && (RXX(t)> RXX(t+1)))
                peak = [peak RXX(t)];
                %location = [location t];
            end
        end
        if(length(peak) > 1)
            MP = max(peak);
            ML = find(RXX == MP);
        end
        M_x=N*(k-1)+ML;
        if(MP >= nguong)
            T0 = ML/Fs;
            F0 = 1/T0;
            if(F0>=75 && F0 <= 450)
                arr(M_x)  = F0;
            end
        end
    end
end
plot(t_x(2:end),arr, '.');
title('Ket Qua');
subplot(2,1,1);
grid;
title('ACF');
xlabel('lags');
ylabel('AutoCorrelation');

function Rxx = autom1(x, k, Fs, K)
    if(k ~= K)
        x1 = x(0.02*Fs*(k-1)+1 : 0.02*Fs*k); % mau khung sau +1 don vi den mau khung sau
    else 
        x1 = x(0.02*Fs*(k-1)+1 : end);
    end
    N = length(x1);
    Rxx = zeros(1, N);
    for m =1: N
        for n=1: N-m
            Rxx(m) = Rxx(m) +x1(n)*x1(n+m);
        end
    end
end




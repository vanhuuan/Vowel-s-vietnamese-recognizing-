[sig,fs] = audioread('lab_female.wav'); %Fetch the signal
%Autocorr_F0 = AutoCorr(sig,fs); %Find F0 in time domain by using AutoCorrelation Function
%Autocorr_F0 = MedSmoothing(Autocorr_F0,7); %Filter the result by with dimension N=7
fft_F0 = mientanso(sig,fs); %Find F0 in frequency domain by using Fast Fourier Transform
fft_F0 = MedSmoothing(fft_F0,7); %Filter the result by Median Smoothing with dimension N=7
figure(1);
subplot(2,1,2)
stem(fft_F0,'filled');
title('Fundamental Frequency F_0(Hz) by Fast Fourier Transform');
ylabel ('F(Hz)');

function [arr_F0] = AutoCorr(inp,fs)
    % The function is used for finding fundamental frequency of a signal in the time domain
    % inp is the input signal
    % fs is sampling frequency
    % The output is arr_F0, which is the array containing calculated F0
    frame_len = round(0.03*fs); % Length of each frame 30ms
    half = round(frame_len/2); % For overlapping frames
    i = 1; % Index of number of elements in the array of F0
    h = hamming(frame_len); % Create a impulse response of the hamming windowd
    disp(h)
    return 
    % Examine every frame of the signal
    for k = 1 : length(inp)/half -1; % (length(inp)/half -1) = number of frames
        width = (k-1)*half + 1:(frame_len + (k-1)*half); % Start index to end index
        frame = h.*inp(width); % Multiplied by the hamming window
        % Use xcorr function to determine the lag of the frame
        [rxx , lag] = xcorr(frame, frame); % Autocorrelation function on the frame
        rxx = rxx(frame_len:end); % Fold the signal
        lag = lag(frame_len:end);
        % Calculate the fundamental frequency F0
        if findpeaks(inp(width),'NPeaks',1,'SortStr','descend') > 0.04 % Remove noise by limiting estimate min
           amplitude = 0.04;
           [y_value,y_peak] = findpeaks(frame,'MinPeakDistance', 285); % Limit min distance between peaks toincrease accuracy
           disp(y_peak);
           lag_a = fs/(y_peak(3) - y_peak(2)); % Find 2 lags by using 3 peaks from the correlation graph
           lag_b = fs/(y_peak(2) - y_peak(1));
            if lag_a<400 && lag_a>80 % Min_f=80Hz and max_f=400Hz
                if lag_b<400 && lag_b>80
                    arr_F0(i) = (lag_a+lag_b)/2;
                    i = i + 1;
                end
            end
        end
    end
end

function [yyy_F0] = mientanso(sig,fs)
    % The function is used for finding fundamental frequency of a signal in spectrum
    % sig is the input signal
    % fs is sampling frequency
    % The output is yyy_F0, which is the array containing calculated F0 in each windown
    N = 32768; %N-point fft
    frame_len = round(0.03*fs); %Length of frame (30ms)
    half = round(frame_len/2); %For overlapping frame
    h = hamming(frame_len); %Hamming window function
    i=1; %Index of element in yyy_F0
    for k = 1 : length(sig)/half -1 %Vong lap cac frame
        range = (k-1)*half + 1:(frame_len + (k-1)*half); %Index of each element in window
        frame = h.*sig(range); %Value of each element in window
        %use FFT function (in Matlab) to analyze the spectrum of the frame
        P2 = abs(fft(frame,N)); %The two-sided spectrum P2
        P1 = P2(1:length(P2)/2+1); %The single-sided spectrum P1
        P1(2:end-1) = 2*P1(2:end-1);
        freq=linspace(0,fs/2,length(P1)); %Spectrum of signal
        if findpeaks(frame,'NPeaks',1,'SortStr','descend') > 0.03 %Remove the noise by limiting the amplitude of signal in frame
            [y_value,y_peak] = findpeaks(P1,freq,'MinPeakHeight',0);%Find peaks (with the certain minimun Height to remove the noise in special cases) in the spectrum of signal
            z_a = y_peak(3) - y_peak(2); %Find 2 F0 by using 3 first peaks in spectrum
            z_b = y_peak(2) - y_peak(1);
            if z<400 && z>80 %Since 80Hz<F0<400Hz, the results will be denied if one of them is outside (80Hz,400Hz)
                if z<400 && z>80
                yyy_F0(i) = (z_a+z_b)/2;%Final result is the mean of 2 F0
                i=i+1; %Put the result into the yyy_F0 and increase the index.
                end
            end
        end
    end
end

function y = MedSmoothing(x,N) %y tin hieu sau khi loc, x la tin hieu truoc khi loc, N la bac can loc
    y = x; %khoi tao gia tri tin hieu y
    temp = 1:N; % tao truc thoi gian tin hieu cua so co chieu dai N
    for k = 1:length(x) % k chay tu 1 den chieu dai cua tin hieu x
        for i = 1:N %i chay tu 1 den N
            if(k < ceil(N/2))
                %xet gia tri tin hieu cua so temp ben trai tin hieu x
                if(i < ceil(N/2)-k+1) % i < phan tu giua cua temp
                    temp(i) = 0;
                else
                    temp(i) = x(k+i-ceil(N/2)); % i > phan tu giua cua temp
                end
            else
                %xet gia tri tin hieu cua so ben phai tin hieu temp
                if(k > length(x) - ceil(N/2) + 1)
                    %i < phan tu giua cua temp
                    if(i < length(x) + ceil(N/2) - k + 1)
                        temp(i) = x(k-ceil(N/2) + i);
                    else
                        temp(i) = 0; %i > phan tu giua cua temp
                    end
                else
                    %xet tin hieu temp chay o giua x.
                    temp(i) = x(k-ceil(N/2)+i);
                end
            end
        end
    %ham median giup tim phan tu trung vi cua day
    y(k) = median(temp);
    end
end
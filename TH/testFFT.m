global Nt Fl Fr fi mat matv slp nfft nds mft mfl nbwl nbws wbwl wbws;

% Parameters
Nt = 0.02; % frame width (in seconds)
fi = 0.01; % frame interval (frame shift in seconds)
mat = 0.04; % MA threshold speech/silence (>=mat: speech, <mat: silence)
matv = 0.13; % MA threshold voiced (>=mat: voiced)
nfft = 1024; % number of fft sampling points
nds = 5; % number of first distances used to evaluate f0
mft = 160; % male/female frequency threshold
mfl = 5; % median filter window length

nbwl = 0.03; % narrowband window length
nbws = 0.01; % narrowband window step
wbwl = 0.005; % wideband window length
wbws = 0.003; % wideband window step

% Constraints
% slp = 0.08;
slp = 0.1; % speech least period (least time range to be considered a speech segment)
           % try to remove virtual speech segments due to noises
Fl = 70; % left boundary of F0
Fr = 400; % right boundary of F0

% Training dataset
[x1,Fs1] = audioread("01MDA.wav");
[x2,Fs2] = audioread("02FVA.wav");
[x3,Fs3] = audioread("03MAB.wav");
[x4,Fs4] = audioread("06FTB.wav");
[x5,Fs5] = audioread("30FTN.wav");
[x6,Fs6] = audioread("42FQT.wav");
[x7,Fs7] = audioread("44MTT.wav");
[x8,Fs8] = audioread("45MDV.wav");

% Training standard results
r1 = [0.45 0.81; 1.53 1.85; 2.69 2.86; 3.78 4.15; 4.84 5.14];
rf1 = [135.5 5.4];
r2 = [0.83 1.37; 2.09 2.60; 3.57 4; 4.76 5.33; 6.18 6.68];
rf2 = [239.7 5.6];
r3 = [1.03 1.42; 2.46 2.8; 4.21 4.52; 6.81 7.14; 8.22 8.5];
rf3 = [115 4.5];
r4 = [1.52 1.92; 3.91 4.35; 6.18 6.6; 8.67 9.14; 10.94 11.33];
rf4 = [202.9 15.5];
r5 = [0.59 0.97; 1.76 2.11; 3.44 3.77; 4.70 5.13; 5.96 6.28];
rf5 = [233.2 11.6];
r6 = [0.46 0.99; 1.56 2.13; 2.51 2.93; 3.79 4.38; 4.77 5.22];
rf6 = [242.7 8.5];
r7 = [0.93 1.42; 2.59 3; 4.71 5.11; 6.26 6.66; 8.04 8.39];
rf7 = [125.7 8.5];
r8 = [0.88 1.34; 2.35 2.82; 3.76 4.13; 5.04 5.5; 6.41 6.79];
rf8 = [177.8 5.7];
rgen = [0 1 0 1 1 1 0 0];

% Train
[t1,f0m1,f0std1,f01,gend1] = plotOutput(x1,Fs1,"01MDA signal",73);
[t2,f0m2,f0std2,f02,gend2] = plotOutput(x2,Fs2,"02FVA signal",73);
[t3,f0m3,f0std3,f03,gend3] = plotOutput(x3,Fs3,"03MAB signal",73);
[t4,f0m4,f0std4,f04,gend4] = plotOutput(x4,Fs4,"06FTB signal",73);
[t5,f0m5,f0std5,f05,gend5] = plotOutput(x5,Fs5,"30FTN signal",73);
[t6,f0m6,f0std6,f06,gend6] = plotOutput(x6,Fs6,"42FQT signal",59);
[t7,f0m7,f0std7,f07,gend7] = plotOutput(x7,Fs7,"44MTT signal",108);
[t8,f0m8,f0std8,f08,gend8] = plotOutput(x8,Fs8,"45MDV signal",108);

% Calculate errors

% Training F0 errors
f0er1 = abs([f0m1 f0std1] - rf1);
f0er2 = abs([f0m2 f0std2] - rf2);
f0er3 = abs([f0m3 f0std3] - rf3);
f0er4 = abs([f0m4 f0std4] - rf4);
f0err1 = f0er1./rf1;
f0err2 = f0er2./rf2;
f0err3 = f0er3./rf3;
f0err4 = f0er4./rf4;
f0err = [mean([f0er1(1) f0er2(1) f0er3(1) f0er4(1)]) mean([f0er1(2) f0er2(2) f0er3(2) f0er4(2)])];

f0er5 = abs([f0m5 f0std5] - rf5);
f0er6 = abs([f0m6 f0std6] - rf6);
f0er7 = abs([f0m7 f0std7] - rf7);
f0er8 = abs([f0m8 f0std8] - rf8);
f0err5 = f0er5./rf5;
f0err6 = f0er6./rf6;
f0err7 = f0er7./rf7;
f0err8 = f0er8./rf8;
f0errt = [mean([f0er5(1) f0er6(1) f0er7(1) f0er8(1)]) mean([f0er5(2) f0er6(2) f0er7(2) f0er8(2)])];

% Training gender errors
genders = [gend1 gend2 gend3 gend4 gend5 gend6 gend7 gend8]; % genders of all signals
correct = sum(abs(genders-rgen) == 0); % number of correct genders
gendaccuracy = correct/length(genders); % gender accuracy rate

% Test & debug
% [f0,y,Fa,n,p] = findF0frame(x1(55*441+1:55*441+882),Fs1);
% [f0,y,Fa,dis] = findF0frame(x2(87*441+1:87*441+882),Fs2);
% [y,f0,dis] = plotF0frame(x2,Fs2,101);

% Calculate and plot result of a signal
% t: ?x2 matrix for speech segments (start ; end) (in seconds)
% f0mean: mean of f0 founded on all frames
% f0std: standard deviation of f0 founded on all frames
% x: signal
% Fs: sampling frequency
function [t,f0mean,f0std,f0,gend] = plotOutput(x,Fs,plottitle,n)
    global mat matv mft nbwl nbws wbwl wbws nfft;
    Tx = 0:(1/Fs):(length(x)-1)*(1/Fs); % time axis
    
    [t,y,Nf] = ma(x,Fs,mat); % find speech segments
    t2 = ma(x,Fs,matv); % find voiced segments
    [f0mean,f0std,f0] = findF0(x,Fs,t2); % find F0
    % evaluate gender
    if f0mean <= mft
        gend = 0;
    else
        gend = 1;
    end
   
    figure
    tiledlayout(4,1)
    nexttile
    % plot input signal
    plot(Tx,x) 
    title(plottitle)
    xlabel("Time (s)")
    ylabel("Magnitude")
    xlim([0 Tx(length(Tx))]);
    
    nexttile
    % plot F0 contour
    plotF0(f0,Tx,t,Fs);
    nexttile
    % plot fft
%     plotF0frame(x,Fs,n);
%     nexttile
    % plot narrowband spectrogram
    spectrogram(x, floor(nbwl*Fs), floor((nbwl-nbws)*Fs), nfft, Fs, 'yaxis');
    ylim([0 8]);
    title("Narrowband spectrogram");
    nexttile
    % plot wideband spectrogram
    spectrogram(x, floor(wbwl*Fs), floor((wbwl-wbws)*Fs), nfft, Fs, 'yaxis');
    ylim([0 8]);
    title("Wideband spectrogram");
end

function [] = plotF0(f0, Tx, t, Fs)
    global Nt fi;
    Ns = floor(Nt/(1/Fs)); % frame width (in samples)
    fis = floor(fi/(1/Fs)); % frame interval (frame shift in samples)
    
    maplot = NaN(1,length(Tx));
    for i=1:length(f0)
        % frame -> time, value is in the middle of frame
        maplot(fis*(i-1)+floor(Ns/2)) = f0(i);
    end
    scatter(Tx, maplot, '.')
    xlim([0 Tx(length(Tx))]);
    ylim([0 400]);
    title("F0")
    ylabel("Frequency (Hz)")
    xlabel("Time (s)")
end

function [y,f0,dis] = plotF0frame(x,Fs,n)
    global Nt fi nds;
    Ns = floor(Nt/(1/Fs)); % frame width (in samples)
    fis = floor(fi/(1/Fs)); % frame interval (frame shift in samples)
    st = (n-1)*fis+1; % start sample index
    ed = (n-1)*fis+Ns; % end sample index
    
    [f0,y,Fa,l,dis] = findF0frame(x(st:ed),Fs);
    [peaks, locs] = findpeaks(y);
    plot(Fa,y)
    xlim([0 Fa(length(Fa))]);
    hold on;
    for i=1:nds+1
        scatter(Fa(locs(i)), peaks(i), 'r')
    end
    hold off;
    title("Magnitude spectrum (Frame " + n + ")");
    xlabel("Frequency (Hz)");
    ylabel("Magnitude");
    legend("FFT", "Selected peaks");
end

function [f0mean, f0std, f0] = findF0(x,Fs,t)
    global Nt fi;
    Ns = floor(Nt/(1/Fs)); % frame width (in samples)
    fis = floor(fi/(1/Fs)); % frame interval (frame shift in samples)
    Nf = floor(length(x)/fis); % Number of frames in signal
    if length(x)-Nf*fis < Ns
        Nf = Nf - ceil(Ns/fis); % Recalculate number of frames - Remove last frames if not enough samples (< frame width)
    end
    
    f0 = NaN(1,Nf);
    % Calculate F0 of frames on founded voiced segments in the signal
    for i=1:size(t,1)
        st = floor(t(i,1)/fi); % start frame index of voiced segment
        ed = floor(t(i,2)/fi); % end frame index of voiced segment
        for j=st:ed
            stf = j*fis+1; % start sample index of frame
            edf = j*fis+Ns; % end sample index of frame
            f0(j) = findF0frame(x(stf:edf),Fs);
        end
    end
    
    % Median filtering 
    f0 = MedianFilter(f0);
    
    f0mean = mean(f0, "omitnan");
    f0std = std(f0, "omitnan");
end

function [f0,y,Fa,locs,dis] = findF0frame(xf,Fs)
    global nfft Fl Fr nds;
    window = hamming(length(xf)); % using hamming window
    y = abs(fft(xf.*window, nfft)); % calculate fft with hamming window
    y = y(1:floor(length(y)/2)); % remove last half of fft due to symmetry values from 0->pi and pi->2pi
    Fa = [0:Fs/nfft:Fs/2-Fs/nfft]; % frequency axis
    
    [peaks, locs] = findpeaks(y);
    % calculate distance between peaks to find F0 if constraints are satisfied
    j = 1;
    for i=1:length(locs)-1
        tmp = Fa(locs(i+1)) - Fa(locs(i));
        if tmp > Fl && tmp < Fr
            dis(j) = tmp;
            j = j + 1;
        end
    end
    
    f0 = mean(dis(1:nds));
end

% Plot MA
% y: array of normalized MA values
% Nf: number of frames in signal
% t: ?x2 matrix for speech segments (start ; end) (in seconds)
function [] = plotMA(y,Nf,t)
    global mat fi;
    plot(1:Nf, y)
    yline(mat,'-.k');
    for i=1:size(t,1)
        % draw vertical lines
        xline(floor(t(i,1)/fi),'r'); % start speech segment
        xline(floor(t(i,2)/fi),'r'); % end speech segment
    end
    title("MA function")
    ylabel("MA value")
    xlabel("Frame")
    legend("Magnitude Average", "Threshold (silence/speech)", "Time boundary (silence/speech)");
    xlim([1 Nf]);
end

% Magnitude Average on signal
% t: ?x2 matrix for speech segments (start ; end) (in seconds)
% y: array of normalized MA values
% Nf: Number of frames in signal
% Ns: frame width (in samples)
% fis: frame interval
% x: signal
% Fs: sampling frequency
% mat: threshold
function [t,y,Nf,fis,Ns] = ma(x,Fs,mat)
    global Nt fi slp;
    Ns = floor(Nt/(1/Fs)); % frame width (in samples)
    fis = floor(fi/(1/Fs)); % frame interval (frame shift in samples)
    Nf = floor(length(x)/fis); % Number of frames in signal
    if length(x)-Nf*fis < Ns
        Nf = Nf - ceil(Ns/fis); % Recalculate number of frames - Remove last frames if not enough samples (< frame width)
    end
    
    y = zeros(1,Nf);
    % Calculate MA on all available frames in the signal
    for i=1:Nf
        y(i) = maf(x,i,Ns,fis);
    end
    % Normalize the result 
    y = Normalize(y);
    
    % Find speech segments
    isSpeech = 0;
    j = 1;
    for i=1:Nf
        % if start speech segment
        if isSpeech == 0 && y(i) >= mat
            isSpeech = 1;
            t(j,1) = (i-1)*fi;
        end
        % if end speech segment
        if isSpeech == 1 && y(i) < mat
            isSpeech = 0;
            t(j,2) = (i-1)*fi;
            % remove speech segments which are too short
            if t(j,2) - t(j,1) <= slp
                t(j,:) = [];
            else
                j = j + 1;
            end
        end
    end
end

% Magnitude Average on frame n
% y: result of MA function on n-th frame
% x: signal
% Ns: frame width (in samples)
% fis: frame interval
function y = maf(x,n,Ns,fis)
    st = (n-1)*fis+1; % start sample index of frame
    ed = (n-1)*fis+Ns; % end sample index of frame
    y = sum(abs(x(st:ed)));
end

function output = MedianFilter(input)
    global mfl;
    hmfl = floor(mfl/2);
    output = input(:);
    for i=hmfl+1:length(input)-hmfl
        output(i) = median(input(i-hmfl:i+hmfl));
    end
end

% Normalize data (point-based)
function [output] = Normalize(input)
    output = input/max(input);
end
global x db
[x,Fs] = audioread("studio_female.wav");
% [x,Fs] = audioread("./TinHieuHuanLuyen/studio_male.wav")
% [x,Fs] = audioread("./TinHieuHuanLuyen/phone_F1.wav")
% [x,Fs] = audioread("./TinHieuHuanLuyen/phone_M1.wav")

fl = 70; % left boundary of F0
fr = 450; % right boundary of F0

% Parameters
Nt = 0.02; % frame width (second)
db = 0.03; % AMDF boundary

global Ns nl nr nfin
Lx = length(x); % number of samples
Ts = 1/Fs; % sample period
Tx = 0:Ts:(Lx-1)*Ts; % point of time for each sample
Ns = Nt/Ts; % number of samples per frame
LN = floor(Lx/Ns); % number of frames
tr = 1/fl; % left boundary of T0
tl = 1/fr; % right boundary of T0
nl = ceil(tl/Ts); % left boundary of lag
nr = floor(tr/Ts); % right boundary of lag
nfin = nl:nr;

% studio_female: x ; first voiced frame: 32 (frame width: 0.02)
% [q,w] = LagAMDF(32)

% AMDF
F = 1:LN;
F0 = zeros(length(F),0);
F0calc = [];
for i=F
    [a,b] = LagAMDF(i);
    F0(i) = 1/(a*Ts);
    % Delete unvoiced and above boundary values
    if ~isnan(F0(i)) && abs(F0(i)-fr)>10 && abs(F0(i)-fl) > 30
        F0calc(length(F0calc)+1) = F0(i);
    else
        F0(i) = NaN;
    end 
end
stem(F,F0,"LineStyle","none",...
    "Marker",".")
F0mean = mean(F0calc)
F0std = std(F0calc)

% Find best lag in frame i
% d: normalized AMDF 
% n: lag value
function [n,d] = LagAMDF(i)
    % Find best lag and its AMDF in frame i
    global Ns nl nr nfin x db
    result = zeros(nr-nl+1,0);
    fst = (i-1)*Ns+1;
    fed = i*Ns;
    input = x(fst:fed);
    min = 1;
    d = NaN; % negative if error 
    n = NaN; % negative if error
    for j=nfin
        result(j-nl+1) = FrameAMDF(input,j) / 100; % normalize by dividing by 100
        if result(j-nl+1) < min
            min = result(j-nl+1);
        else 
            if min < db
                d = min;
                n = j-1;
                break
            end
        end
    end
%     plot(nfin,result)
end

% Short-time AMDF in one frame
% n: integer in [0,N] (n: lag, N: pitch period in samples)
function d = FrameAMDF(xf,n)
    N = length(xf);
    d = sum(abs(xf(1:N-n)-xf(n+1:N)));
end

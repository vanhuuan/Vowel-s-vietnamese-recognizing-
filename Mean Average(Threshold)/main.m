%input_audio = 'lab_female.wav';
input_audio = 'lab_male.wav';
%input_audio = 'studio_female.wav';
%input_audio = 'studio_male.wav';
[y,Fs] = wavread(input_audio);
n_samples = length(y); %so luong mau cua vecto y
duration = n_samples/Fs; % thoi gian thuc cua tin hieu
frame_duration = 0.020; % Moi khung 20ms
min_silence = 0.2;
% Moi file tuong ung mot nguong nhat dinh

    
 
%% Chia khung 100ms
frame_length = Fs * frame_duration; % So mau trong 20ms
n_frames = floor(duration/frame_duration); % So khung
frame_du = mod(n_samples, frame_length); % Khung du
 
frames = reshape( y(1:end-frame_du), frame_length, []);% Matrix frame vs row = so mau trong 20ms va so cot mac dinh = (n_length - frame_du )/so row 
ma = getMA(frames, n_frames, frame_length);
threshold = findThreshold(ma); % tim threshold - duong phan biet khoang lang va giong noi
[A1, silence] = process(frames, threshold,ma, min_silence, n_frames, frame_length, Fs);
processedSignal = A1;
sound(processedSignal, Fs);
%% Ve do thi
draw(processedSignal, y, threshold, ma,silence, frame_du, frame_duration, n_samples, n_frames, Fs);




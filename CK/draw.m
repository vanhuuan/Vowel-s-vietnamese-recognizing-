function [] = draw(processedSignal, y, threshold, ma, silence,frame_du, frame_duration, n_samples, n_frames, Fs)
% Ham ve do thi nhan doi do: Tin hieu ra, tin hieu vao, nguong loc
% Ham ste, Matran silence,so frame du, do dai moi frame, 
% so mau, so frame, Fs
 

subplot(2,1,1);
hold on;
title('Su dung MA de tach am thanh va khoang lang');
xlabel('Thoi gian thuc (s)');
ylabel('Gia tri nang luong da chuan hoa');
 
k = 1:n_frames; % k la so khung
k = k.*(frame_duration); % Chuyen so khung thanh cac moc thoi gian
for i = 1:n_frames
    line(i) = threshold;
end
plot(k, line, 'r','Linewidth',1); % Ve gian do nang luong MA
plot(k, ma, 'g', 'Linewidth',1); % ve nguong di qua gian do MA
 
% Ve cac duong bien
 for i = 2: length(silence)
     j = silence(i)*frame_duration;
     plot([j, j], [0, 1], 'b','LineStyle','-.')
 end 
legend('Nguong','MA', 'Duong bien');
hold off;
 
subplot(2,1,2);
hold on;
xlabel('Thoi gian thuc (s)');
ylabel('Bien do tin hieu');
n = 1:n_samples;
plot(n/Fs, y,'r');
m = 1:n_samples-frame_du;
plot(m/Fs, processedSignal, 'g');
   for i = 1:2:length(silence)-1
        if silence(i) > 0 && silence(i+1) > 0
            disp(i)
            j = silence(i)*frame_duration;
            k = silence(i+1)*frame_duration;
            r = fill([j, k, k, j], [-1, -1, 1, 1], 'r', 'LineStyle','-.') ;
            set(r,'FaceAlpha',0.1) ;
        end
   end
    
% Xu ly cac frame cuoi
if silence(length(silence)) > 0
    j = silence(length(silence))*frame_duration;
    k = n_frames*frame_duration;
    r = fill([j, k, k, j], [-1, -1, 1, 1], 'r','LineStyle','-.') ;
        set(r,'FaceAlpha',0.1);
end
legend('Tin hieu vao','Tin hieu ra','Khoang lang');
hold off;
 
 
end   % end of draw

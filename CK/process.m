function [newframes, silence] = process(frames, threshold, ma, min_silence, n_frames, frame_length, Fs)
% Ham nay nhan doi so la vecto frames va 
% tan so lay mau Fs, khoang thoi gian chia khung,
% muc do loc nhieu
% Tra ve ma tran tin hieu da loc nhieu
% so frame_du, ma, so khung
 
min_silence_samples = min_silence * Fs/frame_length; % tim so mau nho nhat cua doan lang
 
%% Tim giao diem 2 do thi
x = 1: n_frames;  
dif = @(x) (threshold - ma(x)); % tao vecto dif voi gia tri duong la khoang lang, am la tieng noi
 
for i = 2:n_frames
    zx(i) = dif(i)*dif(i-1); 
end  
itsect = find(zx <= 0); % neu tich 2 phan tu lien tiep mang gia tri am thi o do chinh la bien do can tim, itsect chua gia tri la vi tri cac bien do
silence = zeros(1,length(itsect)); % tao vecto silence mau voi so luong = so luong bien do da tim dc
%% Kiem tra am noi hay khoang lang
if length(itsect) > 1 
    for j = 1:length(itsect)-1
        % Nghi van Khoang lang
       if dif(itsect(j)) > 0
           % Neu khoang lang
           if (itsect(j+1) - itsect(j)) > min_silence_samples
               for i = itsect(j):itsect(j+1)
                   silence(j) = itsect(j);
                   silence(j+1) = itsect(j+1);
                   
                    %La khoang lang
                    frames(:,i) = 0;
               end
           end
       end 
    end
    %% Cac frame cuoi
    if dif(itsect(length(itsect))) > 0
        silence(length(itsect)) = itsect(length(itsect));
        for i = itsect(length(itsect)):n_frames
           frames(:,i) = 0;
        end
    end
end
newframes = reshape(frames, [], 1);
end

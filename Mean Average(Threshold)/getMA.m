function [ma] = getMA(frames, n_frames, frame_length)
        % Ham nhan cac doi sola : frames : vector frames,n_frames: do dai cua vector am thanh
        % goc, frame_length: so mau cua 1 frame
        % chieu dai khung
        % Tra ve gian do nang luong MA
        for k = 1:n_frames
            sum = 0;
            for j=1:frame_length
                sum = sum + abs(frames(j,k)); %tong cac phan tu |x[n-m]|
            end
            ma(k) = sum; %MA[k] = tong cac phan tu
        end
        ma_max = max(ma); %tim gia tri lon nhat cua MA
        ma = ma.*(1/ma_max); %chuyen doi bien do cua MA lon nhat thanh 1.
    end % Ket thuc ham getMA


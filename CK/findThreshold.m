function [threshold] = findThreshold(ma) % ham tim threshold cua gian do nang luong MA
        %[a,idx]=unique(ma); % tim overlap cua gian do nang luong MA
        %[ii,ii]=sort(idx); % a(ii) la overlap cua gian do nang luong MA
        f = ma; % vecto f de tim khoang lang
        g = ma; % vecto g de tim giong noi
        Tmin = min(ma);
        Tmax = max(ma);
        T = (Tmin+Tmax)/2;
        i = length(f<T); % so luong gia tri nho hon T
        p = length(g>T); % so luong gia tri lon hon T
        j = -1;
        q = -1;
        while( j~=i || q~=p) % so luong khong doi thi dung vong lap
            sumf = 0; % tong cac phan tu nam trong khoang lang
            sumg = 0; % tong cac phan tu nam trong tieng noi
            for k = 1:length(f)
                    sumf = sumf + max(f(k) - T,0); %tinh tong
                    sumg = sumg + max(T - g(k),0); %tinh tong
            end
            if(sumf > sumg) %neu tong cac phan tu khoang lang > tong cac phan tu tieng noi 
                Tmin = T;        % thi T min = T
            else
                Tmax = T;
            end
            T = (Tmin+Tmax)/2; % Thay doi T
            j = i;
            q = p;
            i = length(f<T); % tinh lai so luong gia tri nho hon T
            p = length(g>T); % tinh lai so luong gia tri lon hon T
        end
        threshold = (T + Tmin)/2;
end
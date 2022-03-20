function [ sample ] = Combine_frame_func( data, frame_length )
    for i = 1 : length(data)
        left = (i - 1)*frame_length/2 + 1;
        right = (i + 1)*frame_length/2;
        sample(left : right) = data(i);
    end
end


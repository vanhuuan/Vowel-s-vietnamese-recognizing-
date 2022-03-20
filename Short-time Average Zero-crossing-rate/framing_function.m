function [frames] = framing_function(data,fs,frame_length)  
    N = length(data);
    num_frames = 2 * floor(N/frame_length) - 1;
    frames(1,:) = data(1 : frame_length);
    for i = 2 : num_frames 
        rightPrev = i*frame_length/2;
        middlePrev = rightPrev - frame_length/2 + 1;
        left = (i - 1)*frame_length/2 + 1;
        right = (i + 1)*frame_length/2;
        halfFramePrev = cat(1, data(middlePrev : rightPrev), zeros(frame_length/2, 1));
        frames(i,:) = data(left : right) + halfFramePrev;
    end
end
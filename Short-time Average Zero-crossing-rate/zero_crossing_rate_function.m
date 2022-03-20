function [ ZCR ] = zero_crossing_rate_function( frames )
    [r, c] = size(frames);
    for i = 1 : r
        frame = frames(i, :);
        ZCR(i) = 0.5*sum(abs(sgn(frame(2:end))-sgn(frame(1:end-1))));
    end
end

    
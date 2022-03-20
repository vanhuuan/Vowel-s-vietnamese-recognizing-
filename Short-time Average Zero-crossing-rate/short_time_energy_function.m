function [ste] = short_time_energy_function(frames)
    [r,c] = size(frames);
    for i = 1 : r
        frame = frames(i, :);
        ste(i) = (1/length(frame))*sum(frame.^2);
    end
end
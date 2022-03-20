function [ result ] = sgn( frame )
    for i=1:length(frame)
        if(frame(i) < 0)
            result(i) = -1;
        else
            result(i) = 1;
        end
    end
end


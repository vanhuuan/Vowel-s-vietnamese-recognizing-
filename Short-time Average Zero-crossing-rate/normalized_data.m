function [ norm_data ] = normalized_data( data )
    norm_data = data/max(abs(data));
    %norm_data = (data - min(data))/(max(data) - min(data));
end


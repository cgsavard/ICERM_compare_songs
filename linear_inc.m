function [y] = linear_inc(b_p_data, params)

% LINEAR_INC is a weight function that weights the y parameter more
% as it increases in the y direction linearly
%
% input: B_P_DATA - nx2 matrix containing the birth and persistence
%                   data
%        PARAMS - contains the min and max value of the persistence
%                 intervals  where you would like to linearly weight
%                 the intervals between
%
% output: Y - weights associated to each interval

min_x = params(1); %x is the interval length
max_x = params(2);

x = b_p_data(:,2); %all d-b values

if x <= min_x
    y = 1;
elseif x >= max_x
    y = 0;
else
    %how they are weighting their points based on how high they are
    y = 1 - (x - min_x)/(max_x - min_x);
end

end
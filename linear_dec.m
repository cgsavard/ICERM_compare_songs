
function [w] = linear_dec(x_y_data, params)

% LINEAR_DEC is a weight function that weights the y parameter less
% as it increases in the y direction linearly
%
% input: X_Y_DATA - nx2 matrix containing the x and y data on the SuPP
%        PARAMS - contains the min and max value of the y data where you 
%        would like to weight the intervals between
%
% output: W - weights associated to each interval

min_u = params(1); %u is the interval length
max_u = params(2);

u = x_y_data(:,2); %all d-b values

if u <= min_u
    w = 1;
elseif u >= max_u
    w = 0;
else
    %how they are weighting their points based on how high they are
    w = 1 - (u - min_u)/(max_u - min_u);
end

end
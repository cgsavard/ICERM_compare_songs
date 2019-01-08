function [y] = gaussian_y(b_p_data, params)

% GAUSSIAN_Y is a weight function that weights the y parameter as a
% gaussian with avg and stdev set by looking at lengths of repeated
% sections on a histogram.
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

%avg and stdev taken from applying gaussian to histo over all lengths
avg = 27.0603;
stdev = 19.0295; 
%apply gaussian in y direction
y = 1/sqrt(2*pi*stdev^2).*exp(-(x-avg).^2./(2*stdev^2));

end
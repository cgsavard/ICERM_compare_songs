function [y] = gaussian_y(x_y_data, params)

% GAUSSIAN_Y is a weight function that weights the y parameter as a
% gaussian with avg and stdev set by looking at lengths of repeated
% sections on a histogram.
%
% input: X_Y_DATA - nx2 matrix containing the x and y data on the SuPP
%        PARAMS - contains the min and max value of the y data where you 
%        would like to weight the intervals between
%
% output: Y - weights associated to each interval

min_u = params(1); %u is the interval length
max_u = params(2);

u = x_y_data(:,2); %all length values

%avg and stdev taken from applying gaussian to histo over all lengths
avg = 27.0603;
stdev = 19.0295; 
%apply gaussian in y direction
y = 1/sqrt(2*pi*stdev^2).*exp(-(u-avg).^2./(2*stdev^2));

end
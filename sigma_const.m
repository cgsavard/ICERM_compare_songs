function [sigx, sigy]= sigma_const(x, y, xnorm)

% SIGMA_CONST uses const. values for the variance at each x and y value.
% The x sigma is constant within each song but dependent on the length
% of the song.
%
% INPUT: X -- the x-values (in SuPP)
%        Y -- the y-values (in SuPP)
%                 
% OUTPUT: SIGX -- the corresponding sigmax given the x-value
%         SIGY -- the corresponding sigmay given the y-value

sigx = 1/xnorm;
sigy = 1;

end
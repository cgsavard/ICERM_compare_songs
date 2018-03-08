function [sigx, sigy] = sigma_parabolic(x, y)

% SIGMA_PARABOLIC uses a parabola for the variance at each x and y value.
%
% INPUT: X -- the x-values
%        Y -- the y-values
%                 
% OUTPUT: SIGX -- the corresponding sigmax given the x-value
%         SIGY -- the corresponding sigmay given the y-value

%Uses a parabola for the sigma values
h = 1/80;
eps = 1/2;
alpha = h*(eps*(1+eps)+1/4);
sigx = -alpha*(x+eps)*(x-(1+eps));
sigy = .5*y+.5;
end
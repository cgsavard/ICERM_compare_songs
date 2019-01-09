function [sigx, sigy] = sigma_parax_liny(x, y, xnorm)

% SIGMA_PARABOLIC uses a parabola for the variance at each x value
% and has a linearly increasing variance in y value.
%
% INPUT: X -- the x-values (in SuPP)
%        Y -- the y-values (in SuPP)
%                 
% OUTPUT: SIGX -- the corresponding sigmax given the x-value
%         SIGY -- the corresponding sigmay given the y-value

%Uses a parabola for the sigma x values
h = 1/80;
eps = 1/2;
alpha = h*(eps*(1+eps)+1/4);
sigx = -alpha*(x+eps)*(x-(1+eps));

%linearly increasing sig y
sigy =  .5*y+.5;

end
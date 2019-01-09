function [mid_norm, norm] = norm_mid(x_y_coords)

% NORM_MID normalizes all start values between 0 and 1 given the max
% and min x
%
% Inputs: X_Y_COORDS - all x,y coords in SuPP
%
% Outputs: X_NORM - new linearly normalized x values
%          NORM - normalization constant for x values

s = x_y_coords(:,1); %start and end
e = x_y_coords(:,2);

mid = (s+e)./2;

mid_max = max(mid);
mid_min = min(mid);

norm = mid_max - mid_min;
mid_norm = (mid-mid_min)./norm;

end



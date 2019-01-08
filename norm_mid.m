function [mid_norm] = norm_mid(b_d_coords)

% NORM_MID normalizes all midpoint values between 0 and 1 given the max
% death and min birth
%
% Inputs: BIRTH - all birth values for a specific song
%
% Outputs: B_LIN_NORM - new linearly normalized birth values

birth = b_d_coords(:,1);
death = b_d_coords(:,2);
d_max = max(death);
b_min = min(birth);

norm = d_max - b_min;
midpoint = (birth + death)/2;

mid_norm = midpoint/norm;

end


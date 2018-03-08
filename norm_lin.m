function [b_lin_norm] = norm_lin(birth)

% NORM_LIN normalizes all birth values linearly between 0 and 1
%
% Inputs: BIRTH - all birth values for a specific song
%
% Outputs: B_LIN_NORM - new linearly normalized birth values

b_max = max(birth);
b_min = min(birth);

b_norm = b_max - b_min;

b_lin_norm = (birth - b_min)/b_norm;

end
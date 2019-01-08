% COLLAPSE_MAPP_XVEC

function [collapse_xvec] = collapse_MaPP_xvec(matrices)

% creates a matrix (collapse_xvec) containing the collapsed x values of 
% the MaPP in each row

num_M = size(matrices,1);
collapse_xvec = [];
for ii = 1:num_M
    collapse_xvec = [collapse_xvec;sum(matrices{ii,1},1)];
end

end
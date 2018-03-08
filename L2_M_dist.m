function [dist_L2] = L2_M_dist(matrices)

% L2_M_DIST computes the L2 (Euclidean) distance between all matrices
% in MATRICES
%
% input: MATRICES - the cell with {n,1} spots filled with matrices
% 
% output: DIST_M - the nxn symmetric L1 distance matrix with the distance
%                  between all pairs of matrices

num_M = size(matrices,1);
D = zeros(num_M, num_M);
for ii = 1:num_M
    for jj = 1:num_M
        A = matrices{ii,1} - matrices{jj,1};
        D(ii,jj) = sqrt(trace(A*A'));
    end 
end

dist_L2 = D;

end
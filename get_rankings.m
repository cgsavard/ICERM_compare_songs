function rank_mat = get_rankings(pwDist, save_flag, save_name)

% GET_RANKINGS converts a pairwise-distance matrix into a ranking
% matrix, which can be used to compute the MAP value for the experiments. 
%
% INPUT: PWDIST -- Matrix of pairwise-distances between songs
%        SAVE_FLAG -- Flag that determines if RANK_MAT will be saved. 
%                     Default is not to save. 
%        SAVE_NAME -- This is the name of the file to be saved. 
%
% OUTPUT: RANK_MAT -- Matrix of rankings by row. 

if nargin < 2
    save_flag = 0;
    save_name = '';
end

% Allocate RANK_MAT
rank_mat = zeros(size(pwDist));

for s = [1:size(pwDist,1)]
    % Take in each row individually and compute the rankings. 
    dist_row = pwDist(s,:);
    
    % Sort distances in descending order
    [~,inds] = sort(dist_row, 'descend');
    
    % Get the rankings for the songs where they are in pwDist by sorting 
    % the indices. 
    [~,rankings] = sort(inds,'ascend');
    
    rank_mat(s,:) = rankings;
end

% Save RANK_MAT if specified. 
if save_flag
    csvwrite(save_name, rank_mat)
end

end
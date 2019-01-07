function [mnn_M] = mutual_nn(dist_M, num_loops)

% MUTUAL_NN looks for mutual nearest neighbors in the data based on the
% distance matrix given.
%
% input: DIST_M - a distance matrix representing the distance between
%                 all songs
%        NUM_LOOPS - indicates how many times you would like to pair
%                    two objects that have not yet been paired the
%                    previous time
% output: MNN_M - mutual nearest neighbor matching vector where the index
%                 is paired with the value at that index

num_songs = size(dist_M,1);
mnn_M = zeros(num_songs,1);

%create new matrix with the diagonal 1000 so that it is not the nn
D = dist_M + eye(num_songs).*10^3;
index_tracker = (1:num_songs)';
%initialize temp variables
new_D = D;
indices = index_tracker;

for loop = 1:num_loops
% set 1 to all indices with min distances
min_M = zeros(size(new_D,1));
for ii = 1:size(min_M,1)
   min_dist = min(new_D(ii,:));
   min_M(ii,:) = (new_D(ii,:) == min_dist).*1;
end

% fill mnn_M with mutual nearest neighbors
compare_M = min_M + min_M';
[row,col] = find(compare_M == 2);
for ii = 1:size(col)   
   col_index = indices(col(ii));
   row_index = indices(row(ii));
   mnn_M(col_index) = row_index;
end

% grab all non-matched indices in mnn_M
indices = index_tracker(find(~mnn_M));

new_D = D(indices,indices);

end

end
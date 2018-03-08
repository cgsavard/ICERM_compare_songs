function [p_matrix] = make_p_diagram(full_matrix_no, full_key)

% MAKE_P_DIAGRAM creates a start-end diagram based on the block lengths
% and start times that are associated with each block present in the
% aligned hierachies
%
% input: FULL_MATRIX_NO - nxm matrix filled with 0 and 1's indicating that
%                         there are blocks of n different interval lengths
%                         and the 1's in the row indicate where the block 
%                         started
%
% output: FULL_KEY - nx1 matrix indicating the length of a block in the 
%                    ith row

%finds row and col associated with 1's
[row,col] = find(full_matrix_no);
num_blocks = length(col);

p_matrix = zeros(num_blocks,2);
for ii = 1:num_blocks
    x = col(ii);
    row_index = row(ii);
    y = x + full_key(row_index);
    p_matrix(ii,:) = [x,y];
end

end
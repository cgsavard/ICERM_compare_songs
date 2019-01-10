function [precision_value, recall_value, recall_wrong] = pr_values(test_M, truth_M)

% PR_VALUES solves for the precision and recall value of the test data
% when compared to the truth data
%
% input: TEST_M - a matrix indicating which songs were matched to which
%                 other songs by inserting a 1 in their shared cell, 
%                 otherwise the cell has a 0
%        TRUTH_M - a matrix containing the truth data in the same format
%                  as above
%
% output: PRECISION_VALUE - out of the songs matched, the percentage
%                           of matches that are correct
%         RECALL_VALUE - out of all the possible matched, the percentage
%                        of matches that are correct

num_correct_matches = length(find(~(test_M - truth_M)));
num_matched = length(find(test_M));
num_songs = length(truth_M);

precision_value = num_correct_matches/num_matched;
recall_value = num_correct_matches/num_songs;

recall_wrong = find((test_M - truth_M))';

end

% COMPARE_SONGS takes in the PIs matrices creates by RUN_SONGS and 
% computes the distance between the matrices and pairs them, then 
% compares the pairings found with the truth data

thresh = '05';
shing = '12';
truth_name = strcat("thresh",thresh,"_shing",shing,"_truth");

% create the truth matrix to compare with
temp = load(strcat("TruthData/",truth_name,".mat"));
truth_M = temp.truth_vec;
%truth_M = make_truth_matrix;

% compute distance matrix from PIs outputted by L1 and L2
%D = L1_M_dist(PIs);
tic
D = L2_M_dist(PIs);
toc

% match these songs using either nearest neighbor or mutual nearest
% neighbor
%matches = nearest_neighbor(D);
%tic
%matches = kmeans(D,4);
matches = mutual_nn(D,1); %int indicates times you would like to pair the unpaired
%toc

% ouput the precision recall values
[p_value, r_value] = pr_values(matches, truth_M)



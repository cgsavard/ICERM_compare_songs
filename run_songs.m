
% RUN_SONGS runs through every song and creates their persistence 
% images which are stored as matrices in the cell songs. COMPARE_SONGS
% should be run next to find distances between the song matrices
% just created

thresh = '05';
shing = '6';
filename = strcat("ScoreData/Thresh",thresh,"_ShingleNumber",shing,"/mazurka*/*.mat");

infiles = dir(char(filename));
num_files = numel(infiles)

songs = cell(0,1);
 
bad_count = 0;
for ii = 1:num_files
    filename = strcat(infiles(ii).folder,"/",infiles(ii).name);
    load (filename);
    [p_matrix] = make_p_diagram(full_matrix_no, full_key);
    if isempty(p_matrix)
        filename
        ' AH empty'
        bad_count = bad_count+1
    else
        songs{end+1,1} = p_matrix;
    end
    
end

tic

res=200; %changes the axes labels
%sig=@sigma_fcn;
sig=@sigma_parabolic;
weight_func=@linear_inc;
params=[0,80]; %[min, max]
norm_fcn = @norm_lin; %Cheb or norm_lin
 %use default setting for hard/soft bounds or specify type=0 or type=1
%[ PIs ] = make_PIs(songs, res, sig, weight_func, params, norm_fcn);

toc

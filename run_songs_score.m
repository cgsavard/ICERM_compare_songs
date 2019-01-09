
% RUN_SONGS runs through every song and creates their persistence 
% images which are stored as matrices in the cell songs. COMPARE_SONGS
% should be run next to find distances between the song matrices
% just created

thresh = '05';
shing = '12';
data = strcat("Thresh",thresh,"_ShingleNumber",shing);

infile = fopen('mazurkas.txt');

num_songs = 52;
songs = cell(0,1);

%truth_vec = 0;
%t_count = 1;

for ii = 1:num_songs
    %good1 = 0;
    %good2 = 0;
    song = string(fgetl(infile));

    filename1 = strcat("ScoreData/",data,"/Expanded/",song);
    load (filename1);
    [p_matrix1] = make_p_diagram(full_matrix_no, full_key);
    if isempty(p_matrix1)
        %good1 = 1;
        filename1
        ' AH empty'
    else
        songs{end+1,1} = p_matrix1;
    end
    
    filename2 = strcat("ScoreData/",data,"/NotExpanded/",song);
    load (filename2);
    [p_matrix2] = make_p_diagram(full_matrix_no, full_key);
    if isempty(p_matrix2)
        %good2 = 1;
        filename2
        ' AH empty'
    else
        songs{end+1,1} = p_matrix2;
    end
    
%     if good1 == 0 && good2 == 0
%         truth_vec(t_count) = t_count+1;
%         truth_vec(t_count+1) = t_count;
%         t_count = t_count+2;
%     elseif good1 == 1
%         good2 = 0;
%     elseif good2 == 1
%         truth_vec(t_count) = 0;
%         t_count = t_count+1;
%     end
end
%truth_vec = truth_vec';

tic

res=200; %changes the axes labels
sig=@sigma_const;
weight_func=@gaussian_y;
params=[0,80]; %[min, max]
norm_fcn = @norm_mid; %norm_mid or norm_lin
type = 1;
 %use default setting for hard/soft bounds or specify type=0 or type=1
[ MaPPs ] = make_MaPPs(songs, res, sig, weight_func, params, norm_fcn, type);

toc
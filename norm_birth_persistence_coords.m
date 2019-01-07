function [ b_p_data, max_b_p, problems] = ...
        norm_birth_persistence_coords(interval_data, norm_fcn )
    
%norm_birth_persistence_coords takes in the interval data as output by the
%duke TDA code (birth-death coordinates) and changes them into
%birth-persistence coordinates with birth normalized from 0 to 1
%
%INPUTS:     interval_data - is a cell array containing the interval data
%            for all of the point clouds in consideration. Each sheet in
%            the cell array corresponds to a different Betti dimension.
%            Each interval set is assumed to be an nX2 matrix with the
%            first column corresponding to the birth times of features and
%            the second column is the death time of each feature. All death
%            times must be greater than the birth time and there must be a
%            finite death time for each feature. 
%            norm_fcn - the function used to normalize the points between
%            0 and 1
%
%OUTPUT:     -b_p_interval_data: This is the modified coordinate data
%            in a cell array. The sheets contain the modified song data.
%            -max_b_p: gives the maximal persistence and maximal
%            birth time across all point clouds for each song. 
%            This information is used to create the boundaries for the
%            persistence images.

[m,n,o]=size(interval_data);
max_persistences=zeros(m,n,o);
max_birth_times=zeros(m,n,o);
birth_persistence=cell(m,n,o);
problems=[];

for k=1:o
for i=1:n
    for j=1:m
        B=interval_data{j,i,k};
        %pulls the song interval data for the (j,i)th point cloud.
        max_persistences(j,i,k)=max(B(:,2)-B(:,1));
        %computes the song persistence (death-birth) for the songs
        max_birth_times(j,i,k)=max(B(:,1));
        %determines that maximal birth time for that songs. We will take 
        %the maximum over all of point clouds to generate non-normalized 
        %song PIs. 
        C=B(:,2)-B(:,1);
        %fcn to normalize births linearly between 0 and 1
        b_norm = norm_fcn(B(:,1));
        birth_persistence{j,i,k}=[b_norm, C];
        %birth-persistence coordinates for song
        D=find(C<0);
        if length(D)>0
            problems=[problems; j,i,k];
        elseif length(D)==0
            problems=problems;
        end

    end
end
song_max_birth(k,1)=1;
%determine the maximum birth time of all song features across the point
%clouds
song_max_persistence(k,1)=max(max(max_persistences(:,:,k)));
%determine the maximum persistence of all song features across the point
%clouds
end
max_b_p=[song_max_birth, song_max_persistence];
b_p_data=birth_persistence;

end
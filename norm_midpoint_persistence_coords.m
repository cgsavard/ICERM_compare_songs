function [ m_p_data, max_coords, problems] = ...
        norm_midpoint_persistence_coords(interval_data, norm_fcn )
    
%norm_midpoint_persistence_coords takes in the interval data as output by the
%duke TDA code (birth-death coordinates) and changes them into
%midpoint-persistence coordinates with birth normalized from 0 to 1
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
%OUTPUT:     -m_p_data: This is the modified coordinate data
%            in a cell array. The sheets contain the modified song data,
%            with the midpoint and persistence.
%            -max_coords: gives the maximal midpoint and maximal
%            persistence across all point clouds for each song. 
%            This information is used to create the boundaries for the
%            persistence images.

[m,n,o]=size(interval_data);
max_persistences=zeros(m,n,o);
max_midpoint_times=zeros(m,n,o);
midpoint_persistence=cell(m,n,o);
problems=[];

for k=1:o
for i=1:n
    for j=1:m
        B=interval_data{j,i,k};
        %pulls the song interval data for the (j,i)th point cloud.
        max_persistences(j,i,k)=max(B(:,2)-B(:,1));
        %computes the song persistence (death-birth) for the songs
        max_midpoint_times(j,i,k)=max(B(:,1));
        %determines that maximal birth time for that songs. We will take 
        %the maximum over all of point clouds to generate non-normalized 
        %song PIs. 
        C=B(:,2)-B(:,1);
        %fcn to normalize midpoint linearly between 0 and 1
        mid_norm = norm_fcn((B(:,2)+B(:,1))/2);
        
        midpoint_persistence{j,i,k}=[mid_norm, C];
        %birth-persistence coordinates for song
        D=find(C<0);
        if length(D)>0
            problems=[problems; j,i,k];
        elseif length(D)==0
            problems=problems;
        end

    end
end
song_max_midpoint(k,1)=1;
%determine the maximum birth time of all song features across the point
%clouds
song_max_persistence(k,1)=max(max(max_persistences(:,:,k)));
%determine the maximum persistence of all song features across the point
%clouds
end
max_coords=[song_max_midpoint, song_max_persistence];
m_p_data=midpoint_persistence;

end
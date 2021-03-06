function [ MaPPs ] = make_MaPPs(interval_data, res, sig, weight_func, params,...
    norm_fcn, type)

%make_MaPPs generates the set of length images for the PH interval data
%stored in the cell array titled interval_data.
%
% INPUTS:    interval__data - is a cell array containing the interval data
%            for all of the point clouds in consideration. Each sheet in
%            the cell array corresponds to a different Betti dimension.
%            Each interval set is assumed to be an nX2 matrix with the
%            first column corresponding to the start times of features and
%            the second column is the end time of each feature. All end
%            times must be greater than the start time and there must be a
%            finite end time for each feature. 
%            res - the desired resolution of the image
%            sig - is the desired variance of the Gaussians used to compute
%            the images
%            weight_func - the name of the weighting function to be
%            called to generate the weightings for the bars. ex
%            params - the set of parameter values needing to be defined for
%            the weighting function. 
%            norm_fnc - the normalization function. Needs to be able to accept the
%            x_y coordinates as an input.
%            type - refers to the declaration of a soft' or 'hard' with
%            respect to the boundaries. type=1 produces hard bounds, type=2
%            produces soft boundaries on the images.
%
%OUTPUTS:    MaPPs - The set of length images generated based on the
%            options specified for the provided interval data.

[ x_y_data, max_x_y, problems ,norm_const] = ...
    norm_length_coords(interval_data, norm_fcn);

%first checks to make sure all the points (start,length) points
%are viable.
if size(problems, 1) > 0
    error('Error: Negative Lengths Present')
elseif size(problems, 1) == 0
    'All positive lengths, continuing'
end

if nargin > 7
    error('Error: Too many input arguments')
elseif nargin == 7
    res = res;
    sig = sig;
    weight_func = weight_func;
    params = params;
    norm_fcn = norm_fcn;
    type = type;
elseif nargin == 6
    res = res;
    sig = sig;
    weight_func = weight_func;
    params = params;
    norm_fcn = norm_fcn;
    type = 1;
elseif nargin == 5
    res = res;
    sig = sig;
    weight_func = weight_func;
    params = params;
    type = 1;
elseif nargin == 4
    error('Error: Incomplete weight function and parameter pair');
elseif nargin == 3 
    res = res;
    sig = sig;
	weight_func = @linear_ramp; %default setting is a linear weighting function
    params = [0, max(max(max_x_y(:,2)))]; %default setting 0 at 0 and 1 at 
    %the maximum length.
    type = 1; %the default boundary setting is hard
elseif nargin == 2
    res = res;
    sig = .5*(max(max(max_x_y(:,2)))/res);%the default setting for the 
    %variance of the gaussians is equal to one half the height of a pixel.
    weight_func = @linear_ramp; %default setting is a linear weighting function
    params = [0, max(max(max_x_y(:,2)))]; %default setting 0 at 0 and 1 at 
    %the maximum length.
    type = 1; %the default boundary setting is hard.
elseif nargin == 1
    res = 25; %default resolution is equal to 25 pixels.
    sig = .5*(max(max(max_x_y(:,2)))/res);%the default setting for the 
    %variance of the gaussians is equal to one half the height of a pixel.
	weight_func = @linear_ramp; %default setting is a linear weighting function
    params = [0, max(max(max_x_y(:,2)))]; %default setting 0 at 0 and 1 at 
    %the maximum length.
    type=1; %the default boundary setting is hard.
end
    
    
if type == 1       
    [ data_images ] = hard_bound_MaPPs( x_y_data, max_x_y, weight_func, ...
        params, res,sig, norm_const);
elseif type == 2
    [ data_images ] = soft_bound_MaPPs( x_y_data, max_x_y, weight_func, ...
        params, res,sig);
end
    MaPPs = data_images;
    
function [ x_y_data, max_x_y, problems, norm_const] = ...
        norm_length_coords(interval_data, norm_fcn )
    
%norm_length_coords takes in the interval data as output by the
%duke TDA code (x-y coordinates) and changes them into
%x-y coordinates with start/midpoint normalized from 0 to 1
%
%INPUTS:     interval_data - is a cell array containing the interval data
%            for all of the point clouds in consideration. Each sheet in
%            the cell array corresponds to a different Betti dimension.
%            Each interval set is assumed to be an nX2 matrix with the
%            first column corresponding to the start times of features and
%            the second column is the end time of each feature. All end
%            times must be greater than the start time and there must be a
%            finite end time for each feature. 
%            norm_fcn - the function used to normalize the points between
%            0 and 1
%
%OUTPUT:     -x_y_data: This is the modified coordinate data
%            in a cell array. The sheets contain the modified song data.
%            -max_x_y: gives the maximal length and maximal
%            start time across all point clouds for each song. 
%            This information is used to create the boundaries for the
%            length images.

[m,n,o] = size(interval_data);
max_y = zeros(m,n,o);
max_x = zeros(m,n,o);
x_y = cell(m,n,o);
problems = [];

for k = 1:o
for i = 1:n
    for j = 1:m
        data = interval_data{j,i,k};
        %pulls the song interval data for the (j,i)th point cloud.
        max_y(j,i,k) = max(data(:,2)-data(:,1));
        %computes the song length (end-start) for the songs
        max_x(j,i,k) = max(data(:,1));
        %determines that maximal start time for that songs. We will take 
        %the maximum over all of point clouds to generate non-normalized 
        %song MaPPs. 
        y = data(:,2)-data(:,1);
        %fcn to normalize starts linearly between 0 and 1
        [x_norm, norm_const] = norm_fcn(data);
        x_y{j,i,k} = [x_norm, y];
        %start-length coordinates for song
        D = find(y<0);
        if length(D) > 0
            problems = [problems; j,i,k];
        elseif length(D) == 0
            problems = problems;
        end

    end
end
%song_max_x(k,1)=1;
song_max_x(k,1) = max(max(max_x(:,:,k)));
%determine the maximum start time of all song features across the point
%clouds
song_max_y(k,1) = max(max(max_y(:,:,k)));
%determine the maximum length of all song features across the point
%clouds
end
max_x_y = [song_max_x, song_max_y];
x_y_data = x_y;

end


function [ data_images ] = hard_bound_MaPPs( x_y_data, max_x_y, ...
        weight_func, params, res, sig, norm_const)
    
%hard_bound_MaPPs generates the MaPPs for a set of point clouds with the
%x_y_data. Hard refers to the fact that we cut the boundaries off hard at
%the maximum values. 
%
%   INPUTS:        x_y_data - start-length points for each of the
%                  point clouds. Each sheet corresponds to a different
%                  Betti dimension.
%                  max_x_y - gives the maximal length and maximal
%                  start time across all point clouds for each song. 
%                  This information is used to create the boundaries for 
%                  the length images.
%                  weight_func - the weighting function the user specified
%                  params - the needed paramteres for the user specified
%                  weight function
%                  res - the resolution (number of pixels). we create
%                  square images with rectangular pixels.
%                  sig - the variance of the gaussians. 
%
%   OUTPUT:        data_images - the set of length images generated
%                  using the selected parameter values. Each sheet
%                  corresponds to the MaPPs generated for different song
%                  interval data.

[m,n,o] = size(x_y_data);
data_images = cell(m,n,o);

for k = 1:o  
song_max_x = 1; %start normalized between 0 to 1
%song_max_y = max_x_y(k,2); 
song_max_y = params(2);

%set up gridding for song
x_stepsize_song = song_max_x/res; %the x-width of a pixel
y_stepsize_song = song_max_y/res; %the y-height of a pixel
grid_values1_song = 0:x_stepsize_song:song_max_x; %must be increasing from zero to max_dist
grid_values2_song = song_max_y:-y_stepsize_song:0; %must be decreasing from max_dist to zero

            for p = 1:m
                for t = 1:n
                song = x_y_data{p,t,k}; %song start length data
                %CHANGES TO THE WIEGHT FUNCTION INPUTS HAPPEN IN THE ROW
                %BELOW
                [weights] = arrayfun(@(row) weight_func(song(row,:), params), 1:size(song,1))';
                [sigmax, sigmay] = arrayfun(@(start, length) sig(start, length, norm_const), song(:,1), song(:,2));
                sigma = [sigmax, sigmay];
                %call the function that makes the image
                [I_song] = grid_SuPP_to_MaPP(song, grid_values1_song, grid_values2_song, sigma,weights);  
                data_images{p,t,k} = I_song;
                end
            end   
end
end


function [ data_images ] = soft_bound_MaPPs( x_y_data, max_x_y, ...
        weight_func, params, res, sig)
    
%soft_bound_MaPPs generates the MaPPs for a set of point clouds with the
%x_y_data. Soft refers to the fact that we add three times the variance
% to the maximal values to determine our boundaries.
%
%   INPUTS:        x_y_data - start-length points for each of the
%                  point clouds. Each sheet corresponds to a different
%                  song
%                  max_x_y - gives the maximal length and maximal
%                  start time across all point clouds for each song. 
%                  This information is used to create the boundaries for 
%                  the length images.
%                  weight_func - the weighting function the user specified
%                  params - the needed paramteres for the user specified
%                  weight function
%                  res - the resolution (number of pixels). we create
%                  square images with rectangular pixels.
%                  sig - the variance of the gaussians. 
%
%   OUTPUT:        data_images - the set of length images generated
%                  using the selected parameter values. Each sheet
%                  corresponds to the MaPPs generated for different song
%                  interval data.

[m,n,o] = size(x_y_data);
data_images = cell(m,n,o);

for k = 1:o    
song_max_x = 1; %start normalized between 0 and 1
song_max_y = max_x_y(k,2);    
    
%set up gridding for song
end_sig = sig(1,0);
x_stepsize_song = (song_max_x+3*end_sig)/res; %the x-width of a pixel
y_stepsize_song = (song_max_y+3*end_sig)/res; %the y-height of a pixel
grid_values1_song = 0:x_stepsize_song:(song_max_x+3*end_sig); %must be increasing from zero to max_dist
grid_values2_song = (song_max_y+3*end_sig):-y_stepsize_song:0; %must be decreasing from max_dist to zero

            for p = 1:m
                for t = 1:n
                song = x_y_data{p,t,k}; %start-length data
                %CHANGES TO THE WIEGHT FUNCTION INPUTS HAPPEN IN THE ROW
                %BELOW
                [weights] = arrayfun(@(row) weight_func(song(row,:), params), 1:size(song,1))';
                %call the funciton that makes the image
                [I_song] = grid_gaussian_bump(song, grid_values1_song, grid_values2_song, sig(t,p),weights);  
                data_images{p,t,k} = I_song;
                end
            end
end
end
%function to linearly interpolate between to get the weighting values for


function [integral_image] = grid_SuPP_to_MaPP(song, grid_values1_song, grid_values2_song, sigma, weights)
    
%grid_SuPP_to_MaPP takes in start-length data for a single point cloud, 
%pixel boundary values,the gaussian variance, and values of the bump
%each gaussian.
%
%Inputs:        BPPairs - is the matrix containing the (start, length)
%               pairs for each interval. This comes from calling the
%               x_y_coordinate funtion.
%               grid_values1_song - is a vector containing the boundaries 
%               for the start values (increasing)
%               grid_values2_song - is a vector containing the boundaries
%               for the start values (decreasing)
%               sigma - is a 1x2 vector with the x and y standard
%               deviations.
%               weights - vector containing the weighting value for each
%               interval as determined by the user specified weighting
%               function and parameters.
%
%Outputs:       integral_image - is the image computed by discreting based
%               on the values contained in grid_values and summing over all
%               thedifferent gaussians centered at each point in the
%               start-length interval data.

max_bar_length = grid_values2_song(1);

%keyboard
[X,Y] = meshgrid(grid_values1_song,grid_values2_song);
XX = reshape(X,[],1);
YY = reshape(Y,[],1);
for k = 1:size(song,1) %looks at each point in the song
    M = weights(k);
    s = sigma(k,:);
    AA = (M)*mvncdf([XX YY], song(k,:), s);
    AA = reshape(AA,[],length(grid_values1_song));
    ZZ(:,:,k) = (-AA(2:end,2:end)-AA(1:end-1,1:end-1)+AA(2:end, 1:end-1)+AA(1:end-1,2:end));
    % The above line is implementing the procedure explained in https://en.wikipedia.org/wiki/Summed_area_table
end

%finding the max valume under the gaussian per pixel
integral_image = max(ZZ,[],3);

end


end

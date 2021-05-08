%% Introduction
% Computes the average distance (location and orientation) of source states
% used in the training and test sets.

%% Setup
addpath(genpath('./functions'))
clear
clc
pre_process_figure();

%% Compute distances
files = {
    './config/locations/train_sources_res0.09.mat', 
    './config/locations/train_sources_res0.05.mat',
    './config/locations/train_sources_res0.03.mat',
%     './config/locations/train_sources_res0.01.mat'
};
min_distances = [0.09, 0.05, 0.03, 0.01];
avg_min_loc = zeros(size(min_distances));
avg_min_azi = zeros(size(min_distances));

idx = 1
for f = files'   
    load(f{1});
    sources = train_sources;
    clear train_sources
    
    avg_min_loc(idx) = location_distances(sources(:, 1:2));
    avg_min_azi(idx) = orientation_distances(sources(:, 3));
    idx = idx + 1;
end

res = table(min_distances', avg_min_loc', avg_min_azi',...
            'VariableNames', {'sample_distance', 'location', 'orientation'});
disp(res)

function avg_min_d = location_distances(X)

d = pdist(X, 'euclidean'); % condensed
min_d = mink(squareform(d), 2, 2);
avg_min_d = mean(min_d(:, 2));

end


function avg_min_d = orientation_distances(X)

d = pdist(X, @compute_angle_difference); % condensed
min_d = mink(squareform(abs(d)), 2, 2);
avg_min_d = mean(min_d(:, 2));

end

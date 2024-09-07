% Script to calculate footprint areas, infer repo volumes, and estimate
% fluxes
clear
close all

DataPath = '../dat/20220214/';
FileName = 'sbapts_02142022putm.csv';
OutDir = '../out/20220214/';
mkdir(OutDir)

% load table
refTable = readtable([DataPath FileName]);
keyTable = readtable([OutDir 'key_table.csv']);

load([OutDir 'height_fit_struct.mat'],'height_fit_struct')

% specify number of bootstraps
nBoots = 100;
angle_ind = 2;

% initialize data structure
vol_capacity_array = NaN(size(keyTable,1),nBoots);
sa_capacity_vec = NaN(size(keyTable,1),1);
mean_x_vec = NaN(size(keyTable,1),1);
mean_y_vec = NaN(size(keyTable,1),1);

% iterate through repositories
% wb = waitbar(0, 'Estimating repository volume');
tic
parfor r = 1:size(keyTable,1)  
%     waitbar(r/size(keyTable,1), wb);
    
    repo_footprint_filter = keyTable.rep_id(r) == refTable.ORIG_FID;           

    % extract footprint and calculate total volume    
    xfp = refTable.POINT_X(repo_footprint_filter);
    yfp = refTable.POINT_Y(repo_footprint_filter);
                 
    % calculate mean position
    mean_x_vec(r) = nanmean(xfp);
    mean_y_vec(r) = nanmean(yfp);
    
    [height_array, dist_array, repo_mask, x_bounds, y_bounds,max_height_vals] = ...
                      vol_calculations(xfp,yfp,height_fit_struct(angle_ind),[],[],[],nBoots,[]);

    vol_capacity_array(r, :) = reshape(sum(sum(height_array,1),2),1,[]);
    sa_capacity_vec(r) = sum(repo_mask(:));
        
end    
% delete(wb)
toc

% generate summary table
volume_table = struct;
volume_table.names_clean = keyTable.names_clean;
volume_table.names_raw = keyTable.names_raw;
volume_table.region_id = keyTable.region_id;
volume_table.repo_id = keyTable.rep_id;
volume_table.date_num = keyTable.date_num;
volume_table.vol_capacity_mean = nanmean(vol_capacity_array,2);
volume_table.vol_capacity_ste = nanstd(vol_capacity_array,[],2);
volume_table.sa_capacity = sa_capacity_vec;
volume_table.mean_x_pos = mean_x_vec(:,1);
volume_table.mean_y_pos = mean_y_vec(:,1);

volume_table = struct2table(volume_table);

writetable(volume_table,[OutDir 'repo_volume_estimates.csv'])
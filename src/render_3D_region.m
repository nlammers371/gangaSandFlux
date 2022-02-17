% Script to make visualization of cluster of repos
clear
close all

addpath(genpath('utilities'));
RawDataPath = '../dat/';
ProcessedDataPath = '../out/';
mkdir(ProcessedDataPath);
FigDir = '../fig/';
mkdir(FigDir)

RepFootprintFileCell = {'all_gangareps_4_23.csv','GangaReps_0910pts.csv'};q
suffix_cell = {'2020','2010'};
RepoIDCell = {887:920,153:184};
angle_rep_index = 2;

for i = 1:length(RepoIDCell)
    
    RepFootprintFile = RepFootprintFileCell{i};


    % specify IDs for region to visualize
    regions_to_plot = RepoIDCell{i};
    
    % load table
    footprintTable = readtable([RawDataPath RepFootprintFile]);

    % load reference set for repo heights
    load([ProcessedDataPath 'height_fit_struct.mat'],'height_fit_struct')
    angle_rep_vec = [height_fit_struct.angle_rep];
    angle_of_repose = angle_rep_vec(angle_rep_index);

    % extract vectors
    x_vec = footprintTable.POINT_X;
    y_vec = footprintTable.POINT_Y;

    repo_id_vec = round(footprintTable.ORIG_FID);
    repo_filter = ismember(repo_id_vec,regions_to_plot);

    % apply filter
    repo_id_vec_ft = repo_id_vec(repo_filter);
    x_vec_ft = x_vec(repo_filter);
    y_vec_ft = y_vec(repo_filter);

    % renormalize position data
    norm_factor = 1e1;

    min_x = floor(min(x_vec_ft)/norm_factor)*norm_factor;
    min_y = floor(min(y_vec_ft)/norm_factor)*norm_factor;

    max_x = ceil((max(x_vec_ft)-min_x)/norm_factor)*norm_factor;
    max_y = ceil((max(y_vec_ft)-min_y)/norm_factor)*norm_factor;

    x_shifted = (x_vec_ft - min_x);
    y_shifted = (y_vec_ft - min_y);

    % initialize empty frame to populate
    region_fp_frame = zeros(max_y,max_x);
    region_dist_frame = zeros(max_y,max_x);
    region_height_frame = zeros(max_y,max_x);

    rng(102);

    for r = 1:length(regions_to_plot)

        repo_footprint_filter = repo_id_vec_ft == regions_to_plot(r);           

        % extract footprint and calculate total volume    
        xfp = x_shifted(repo_footprint_filter);
        yfp = y_shifted(repo_footprint_filter);

        % generate mast
        repo_mask = poly2mask(xfp, yfp, max_y, max_x);    
        region_fp_frame = region_fp_frame | repo_mask;

        % calculate SA
        SA = sum(repo_mask(:));

        % generate distance map
        dist_array = bwdist(~repo_mask);
        region_dist_frame(repo_mask==1) = dist_array(repo_mask==1);

        % calculate max height (assume deterministic for now) 
        mh_mean = height_fit_struct(angle_rep_index).x_fit(1) * SA ./ (height_fit_struct(angle_rep_index).x_fit(2) + SA);
        mh_sigma = height_fit_struct(angle_rep_index).mdl_se;
        min_height = min(height_fit_struct(angle_rep_index).max_height_vec);
        height_dist = makedist('normal', mh_mean, mh_sigma);
        height_dist = truncate(height_dist,min_height,Inf);    
        height_val = random(height_dist,1,1);

        % calculate repo capacity   
        height_array = zeros(size(dist_array,1),size(dist_array,2));    
        dist_to_top = height_val/tand(angle_of_repose);        
        height_array(dist_array>=dist_to_top) = height_val;
        height_array(dist_array<dist_to_top) = dist_array(dist_array<dist_to_top) ...
                                                        * tand(angle_of_repose);

        region_height_frame(repo_mask==1) = height_array(repo_mask==1);

    end

    %%% Make figure
%     close all
    surf_fig = figure;
    cmap = colormap('pink');
    cmap = cmap(100:end,:);
    colormap(cmap);
    height_slice = imgaussfilt(region_height_frame,2); % apply mild smoothing
    surf(height_slice,'EdgeAlpha',0.1)        
    % xlabel('x distance (meters)')
    % ylabel('y distance (meters)')
    zlabel('inferred heigh (meters)')

    zlim([0 30])
    xlim([0 size(height_slice,2)])
    ylim([0 size(height_slice,1)])
    if i == 1
        view(-45,45)
    else
        view(-45,45)
    end
    set(gca,'Fontsize',10)

    saveas(surf_fig,[FigDir 'ROI_reconstruction_' suffix_cell{i} '.png'])
    saveas(surf_fig,[FigDir 'ROI_reconstruction_' suffix_cell{i} '.pdf'])
end



clear
close all

RawDataPath = '../dat/';
ProcessedDataPath = '../out/';
mkdir(ProcessedDataPath);
RepFootprintFiles = {'all_gangareps_4_23.csv','GangaReps_2015pts.csv','GangaReps_0910pts.csv'};
FigDir = '../fig/';
mkdir(FigDir)

master_capacity_struct = struct;
n_boots = 100;

for f = 1:length(RepFootprintFiles)

    % load table
    footprintTable = readtable([RawDataPath RepFootprintFiles{f}]);

    % load reference set for repo heights
    load([ProcessedDataPath 'height_fit_struct.mat'],'height_fit_struct')
    angle_rep_vec = [height_fit_struct.angle_rep];

    % extract vectors
    x_vec = footprintTable.POINT_X;
    y_vec = footprintTable.POINT_Y;
    
%     area_vec = footprintTable.area;
    repo_id_vec = footprintTable.ORIG_FID;
    repo_id_index = unique(repo_id_vec);
    % renormalize position data
    norm_factor = 1e4;

    min_x = min(x_vec);
    min_y = min(y_vec);

    max_x = ceil((max(x_vec)-min_x)/norm_factor)*norm_factor;
    max_y = ceil((max(y_vec)-min_y)/norm_factor)*norm_factor;

    x_shifted = (x_vec - min_x);
    y_shifted = (y_vec - min_y);


    vol_capacity_vec = NaN(n_boots,length(repo_id_index),length(angle_rep_vec));    
    sa_fp_vec = NaN(length(angle_rep_vec),length(repo_id_index));

    % struct for footprints and predicted pile capacity
    repo_capacity_struct = struct;
    % iterate through each repository 
    for r = 1:length(repo_id_index)

        repo_footprint_filter = repo_id_vec == repo_id_index(r);           

        % extract footprint and calculate total volume    
        xfp = x_shifted(repo_footprint_filter);
        yfp = y_shifted(repo_footprint_filter);

        for a = 1:length(angle_rep_vec)                            

           [height_array, dist_array, repo_mask, x_bounds, y_bounds,max_height_vals] = ...
                          vol_calculations(xfp,yfp,height_fit_struct(a),[],[],[],n_boots,[]);

            vol_capacity_vec(:,r,a) = reshape(sum(sum(height_array,1),2),[],1);
            sa_fp_vec(a,r) = sum(repo_mask(:));

            % record repo info               
            repo_capacity_struct(r).example_height_array(:,:,a) = height_array(:,:,1); % just take a random example


        end    
        repo_capacity_struct(r).repo_dist_array = dist_array;
        repo_capacity_struct(r).sa_array = repo_mask; 
        repo_capacity_struct(r).x_points = xfp;
        repo_capacity_struct(r).y_points = yfp;
        repo_capacity_struct(r).id = r;    

    end
    master_capacity_struct(f).repo_capacity_struct = repo_capacity_struct;
    master_capacity_struct(f).vol_capacity_array = vol_capacity_vec;
    master_capacity_struct(f).sa_fp_array = sa_fp_vec;
    master_capacity_struct(f).repo_id_index = repo_id_index;    
    master_capacity_struct(f).n_boots = n_boots;
    master_capacity_struct(f).source_file = RepFootprintFiles{f};
    master_capacity_struct(f).source_data = footprintTable;
end    
    
save([ProcessedDataPath 'master_capacity_struct.mat'],'master_capacity_struct')
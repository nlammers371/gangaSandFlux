clear
close all

DataPath = '../dat/';
OutPath = '../out/';
mkdir(OutPath);
FileName = 'mirz_1516flux.csv';
RefName = 'Rep_Stats_all.xlsx';
FigDir = '../fig/';
mkdir(FigDir)

% load table
pointTable = readtable([DataPath FileName]);
infoTable = readtable([DataPath RefName]);

%% calculate height-area function for height sampling
area_ref_vec = infoTable.base_area;
height_fun = @(x) x(1) * area_ref_vec ./ (x(2) + area_ref_vec);
max_obs_height = max(pointTable.area);
% make filter for points that are in the relevant range
rel_area_filter = area_ref_vec <=max_obs_height;
% make new height field
infoTable.height_35 = infoTable.height_30*tand(35)/tand(30);
height_fit_struct = struct;
angle_rep_vec = [25 30 35];
area_index = linspace(min(area_ref_vec),2*max(area_ref_vec));

for a = 1:length(angle_rep_vec)
    eval(['ob_fun = @(x) infoTable.height_' num2str(angle_rep_vec(a)) '- height_fun(x);'])
    [x_fit,~,resid] = lsqnonlin(ob_fun,[10 1e4],[0 0],[50,Inf]);    
    height_fit_struct(a).x_fit = x_fit;
    height_fit_struct(a).angle_rep = angle_rep_vec(a);
    height_fit_struct(a).mdl_se = sqrt(sum(resid(rel_area_filter).^2)/sum(rel_area_filter));
    height_fit_struct(a).max_height_vec = x_fit(1).*area_index./(x_fit(2) + area_index);
end    
    
save([OutPath 'height_fit_struct.mat'],'height_fit_struct')

% extract vectors
x_vec = pointTable.POINT_X;
y_vec = pointTable.POINT_Y;
area_vec = pointTable.area;
master_id_vec = pointTable.ORIG_FID;

% renormalize position data
norm_factor = 1e4;

min_x = min(x_vec);
min_y = min(y_vec);

max_x = ceil((max(x_vec)-min_x)/norm_factor)*norm_factor;
max_y = ceil((max(y_vec)-min_y)/norm_factor)*norm_factor;

snapshot_index = unique(master_id_vec);

x_shifted = (x_vec - min_x);
y_shifted = (y_vec - min_y);


% let's do some basic string parsing to extract the dateinfo
name_vec = pointTable.name;
region_string_cell = cell(size(snapshot_index));
date_vec = NaN(size(name_vec));
date_cell = cell(size(snapshot_index));
footprint_flag_vec = false(size(name_vec));
for n = 1:length(snapshot_index)
    snap_filter = master_id_vec == snapshot_index(n);
    first_ind = find(snap_filter,1);
    string_temp = name_vec{first_ind};
    dash_ind = strfind(string_temp,'_');
    region_string_cell{n} = string_temp(1:dash_ind-1);
    date_string = string_temp(dash_ind+1:end);
    if strcmpi(date_string,'perimeter')
        footprint_flag_vec(snap_filter) = true;
        date_cell{n} = NaN;        
    else
        Y = str2double(['20' date_string(1:2)]);
        M = str2double(date_string(3:4));
        D = str2double(date_string(5:6));
        dt = datetime(Y,M,D);
        date_cell{n} = dt;
        date_vec(snap_filter) = datenum(dt);
    end
end

% generate region id vec
unique_regions = unique(region_string_cell);
repo_id_index = 1:length(unique_regions);
repo_id_vec = NaN(size(master_id_vec));
for n = 1:length(snapshot_index)
    snap_filter = master_id_vec == snapshot_index(n);
    % find region index
    region_ind = find(strcmp(region_string_cell{n},unique_regions));
    repo_id_vec(snap_filter) = region_ind;
end

% let's generate vol estimates and footprints for each slice/time point

% assume a single angle of repose for now
angle_of_repose = 25;

% initialize array to store results
date_index = unique(date_vec(~isnan(date_vec)));

% struct for footprints and predicted pile
n_boots = 100;
master_flux_struct = struct;

for a = 1:length(angle_rep_vec)  
    repo_flux_struct = struct;
    
    vol_flux_array = NaN(length(date_index),length(repo_id_index),n_boots);
    vol_fp_vec = NaN(n_boots,length(repo_id_index));
    sa_flux_array = NaN(length(date_index),length(repo_id_index));
    sa_fp_vec = NaN(1,length(repo_id_index));

    % iterate through each repository 
    for r = 1:length(repo_id_index)

        repo_footprint_filter = footprint_flag_vec == 1 & repo_id_vec == repo_id_index(r);    

        % calculate the expected maximum height
%         sa_input = mean(area_vec(repo_footprint_filter));

        % extract footprint and calculate total volume    
        xfp = x_shifted(repo_footprint_filter);
        yfp = y_shifted(repo_footprint_filter);

       [height_array, dist_array, repo_mask, x_bounds, y_bounds,max_height_vals] = ...
                      vol_calculations(xfp,yfp,height_fit_struct(a),[],[],[],n_boots,[]);
%         (x_vec,y_vec,fit_struct,x_bounds,y_bounds,SA,n_boots, max_height_vals)    
        vol_fp_vec(:,r) = reshape(sum(sum(height_array,1),2),[],1);
        sa_fp_vec(r) = sum(repo_mask(:));
        
        % record repo info
        repo_flux_struct(r).sa_array = repo_mask;
        repo_flux_struct(r).repo_dist_array = dist_array;
        repo_flux_struct(r).repo_height_array = height_array(:,:,1);
        repo_flux_struct(r).x_points = xfp;
        repo_flux_struct(r).y_points = yfp;

        % now iterate through dates and calculate footprint for each
        repo_flux_struct(r).sa_stack = NaN(size(height_array,1),size(height_array,2),length(date_index));
        repo_flux_struct(r).dist_stack = NaN(size(height_array,1),size(height_array,2),length(date_index));
        repo_flux_struct(r).height_stack = NaN(size(height_array,1),size(height_array,2),length(date_index));
        repo_flux_struct(r).id = r;

        for d = 1:length(date_index)

            % filter
            dt_filter = date_index(d)==date_vec & repo_id_vec == repo_id_index(r);  

            if sum(dt_filter)>1

                xdt = x_shifted(dt_filter);
                ydt = y_shifted(dt_filter);

                % calculate volume
                [height_array_dt, dist_array_dt, repo_mask_dt] = ...
                      vol_calculations(xdt,ydt,height_fit_struct(a),x_bounds,y_bounds,sa_fp_vec(r),n_boots,max_height_vals); 

                % record
                vol_flux_array(d,r,:) = sum(sum(height_array_dt,1),2);
                sa_flux_array(d,r) = sum(repo_mask_dt(:));
                a_rep = mean(area_vec(dt_filter));

                % add to struct
                repo_flux_struct(r).sa_stack(:,:,d) = repo_mask_dt;
                repo_flux_struct(r).dist_stack(:,:,d) = dist_array_dt;
                repo_flux_struct(r).height_stack(:,:,d) = height_array_dt(:,:,1);
            end
        end 
    end    

    % fill in the gaps assuming vol = 0 for any missing observations between
    % first and last detection
    index_vec = (1:size(vol_flux_array,1))';
    iter = 1;
    for r = 1:size(vol_flux_array,2)
        vol_vec = vol_flux_array(:,r,1);
        first_i = find(~isnan(vol_vec),1);
        last_i = find(~isnan(vol_vec),1,'last');
        nan_flags = index_vec>=first_i&index_vec<=last_i& isnan(vol_vec);

        iter = iter + sum(nan_flags);

        vol_flux_array(nan_flags,r,:) = 0;
        sa_flux_array(nan_flags,r) = 0;
    end    
    
    master_flux_struct(a).repo_flux_struct = repo_flux_struct;
    master_flux_struct(a).vol_flux_array = vol_flux_array;
    master_flux_struct(a).sa_flux_array = sa_flux_array;
    master_flux_struct(a).vol_fp_array = vol_fp_vec;
    master_flux_struct(a).sa_fp_array = sa_fp_vec;
    master_flux_struct(a).repo_id_index = repo_id_index;
    master_flux_struct(a).date_index = date_index;
    master_flux_struct(a).n_boots = n_boots;
end

save([OutPath 'master_flux_struct.mat'],'master_flux_struct')
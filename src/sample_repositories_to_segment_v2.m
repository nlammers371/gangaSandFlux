% Script to randomly select subset of repositories to segment in detail for
% volume flux analysis
clear
close all

DataPath = '../dat/20220214/';
FileName = 'sbapts_02142022putm.csv';
OutDir = '../out/20220214/';
mkdir(OutDir)

% load table
refTable = readtable([DataPath FileName]);

% First, clean up location names
refTableClean = refTable;
[NameIndexRaw, ia, ic] = unique(refTable.Name);

x_vec = refTable.POINT_X;
y_vec = refTable.POINT_Y;

NameIndexClean = NameIndexRaw;

% calculate average x and y positions
x_mean_vec = grpstats(x_vec,ic,'mean');
y_mean_vec = grpstats(y_vec,ic,'mean');

tic
% clean up name formats
for i = 1:length(NameIndexRaw)
    name = NameIndexRaw{i};
    
    % account for random bullshit
    if contains(name,'_E_')
        name = strrep(name,'_E_','_');
        name = ['e' name];        
    elseif contains(name,'_C_')
        name = strrep(name,'_C_','_');
        name = ['c' name];        
    elseif contains(name,'_W_')
        name = strrep(name,'_W_','_');
        name = ['w' name];        
    elseif contains(name,'_S_')
        name = strrep(name,'_S_','_');
        name = ['s' name];        
    elseif contains(name,'_N_')
        name = strrep(name,'_N_','_');
        name = ['n' name];  
    elseif contains(name,'__')
        name = strrep(name,'__','_');
    end
    
    underscores = strfind(name,'_');
    % compress cardinal info
    if underscores(1)==2
      name = [name(1) name(3:end)];
    end
    
    % append boundary info if not there
    underscores = strfind(name,'_');
    if length(underscores)==2
        name = [name '_n'];
    end
    NameIndexClean{i} = name;
end
toc
% initialize new structure
[NameIndexClean,map_to] = unique(NameIndexClean);
NameIndexRaw1 = NameIndexRaw(map_to);

RepoTable.RepoNames = NameIndexClean;

% now parse the names
regionCell = cell(size(NameIndexClean));
repoCell = cell(size(NameIndexClean));
dateStringCell = cell(size(NameIndexClean));
dateNumVec = NaN(size(NameIndexClean));
subIDVec = NaN(size(NameIndexClean));
uniqueIDVec = (1:length(NameIndexClean))';
tic
for i = 1:length(NameIndexClean)
    name = NameIndexClean{i};
    % find underscores   
    underscores = strfind(name,'_');
    % get region name
    regionCell{i} = name(1:underscores(1)-1);
    % extract date
    dateStringCell{i} = name(underscores(1)+1:underscores(2)-1);
%     dt = datetime(date_string,'InputFormat','MMddyyyy');
%     dateVec(i) = dt;
    % repo sub ID
    subIDVec(i) = str2double(name(underscores(2)+1:underscores(3)-1));
    % add concatenated repo ID field
    repoCell{i} = [name(1:underscores(1)-1) name(underscores(2)+1:underscores(3)-1)];
end    


% get list of unique dates
[unique_dates, map_to, map_from] = unique(dateStringCell);
for u = 1:length(unique_dates)
    dt = datenum(datetime(unique_dates{u},'InputFormat','MMddyyyy'));
    dateNumVec(map_from==u) = dt;
end    
toc
% Now generate vector of region IDs and a master ID vector
% dateNumVec = datenum(dateVec);
[region_id_index, ~, region_id_vec] = unique(regionCell);

% now get array with unique region-date combinations
[region_date_array,~,rd_map_vec] = unique([region_id_vec dateNumVec],'rows');

% look at consecutive differences
region_date_diffs = diff(region_date_array,1,1);
region_date_diffs(region_date_diffs(:,1)~=0,2) = Inf;

% calculate average x and y positions
x_region_mean_vec = grpstats(x_mean_vec,rd_map_vec,'mean');
y_region_mean_vec = grpstats(y_mean_vec,rd_map_vec,'mean');
[~,region_pca_array] = pca([x_mean_vec y_mean_vec]);

% Identify duples with sufficiently close time points
max_dt = 15; % maximum temporal separation
n_repos_per_region = 15;
nBins = 5;
total_samples = 1e3;

% apply dt filter and generate sampling weights
dt_filter_vec = find([region_date_diffs(:,2)<=max_dt ; false]);
x_region_mean_vec_ft = x_region_mean_vec(dt_filter_vec);
y_region_mean_vec_ft = y_region_mean_vec(dt_filter_vec);
date_vec_ft = region_date_array(dt_filter_vec,2);
region_id_vec_ft = region_date_array(dt_filter_vec,1);
dt_vec_filtered = region_date_diffs(dt_filter_vec,2);

% do pca
[~, dist_pca_array] = pca([x_region_mean_vec_ft y_region_mean_vec_ft]);

% get counts
[count_array,~,~,xID,yID] = histcounts2(date_vec_ft,dist_pca_array(:,1),nBins);
linIndex = sub2ind(size(count_array),xID,yID);

% generate weight vector 
sample_weight_vec = (1./dt_vec_filtered).^2;%.* 1./count_array(linIndex);

%% perform actual sampling
sampled_flags = false(size(sample_weight_vec));
sample_counter = 1;
iter_counter = 1;
sample_struct = struct;
rng(147);

while ~all(sampled_flags) && sample_counter < total_samples
  
    % randomly draw region/date id
    samp_options = find(~sampled_flags);
    if length(samp_options) > 1
        samp_id = randsample(samp_options,1,true,sample_weight_vec(samp_options));
    else
        samp_id = samp_options;
    end
    
    % update tracking vector
    sampled_flags(samp_id) = true;
    
    % identify candidate repositories
    region_id = region_date_array(dt_filter_vec(samp_id),1);
    date_id = region_date_array(dt_filter_vec(samp_id),2);
    date_id2 = region_date_array(dt_filter_vec(samp_id)+1,2);
    repo_ids = find(region_id_vec == region_id & ...
                       dateNumVec == date_id);
                  
    % randomly select repos to use
    repo_rank_list = randsample(1:length(repo_ids),min([length(repo_ids) 2*n_repos_per_region]),false);
    repo_id_list = repo_ids(repo_rank_list);
    
    for r = 1:length(repo_rank_list)
      
        repo_id = repo_id_list(r);
        
        % add info
        sample_struct(iter_counter).region_name = region_id_index{region_id};
        sample_struct(iter_counter).region_id = region_id;
        sample_struct(iter_counter).repo_rank_vec = r;
        sample_struct(iter_counter).primary_sample_flag = 1*(r < n_repos_per_region); 
        sample_struct(iter_counter).date_1 = dateStringCell{repo_id};
        sample_struct(iter_counter).date_2 = datestr(date_id2,'mmddyyyy');
        sample_struct(iter_counter).repo_name_raw = NameIndexRaw1(repo_id);
        sample_struct(iter_counter).repo_name_clean = NameIndexClean(repo_id);
        sample_struct(iter_counter).mean_x_pos = x_mean_vec(repo_id);
        sample_struct(iter_counter).mean_y_pos = x_mean_vec(repo_id);
        sample_struct(iter_counter).mean_pca_pos = region_pca_array(repo_id);                
        sample_struct(iter_counter).date_1_num = dateNumVec(repo_id);

        iter_counter = iter_counter + 1;
        sample_counter = sample_counter + 1*(r < n_repos_per_region);
    end
end  

% Convert to table
sample_table = struct2table(sample_struct);
% reorder by region-date-rank order
sample_table = sortrows(sample_table,[2 12 3]);

writetable(sample_table,[OutDir 'repo_sample_list.csv'])
%% Look at distribution of selected regions
samp_pos_vec = [sample_struct.mean_pca_pos]/1e4;
samp_date_vec = [sample_struct.date_1_num]/1e5;
samp_prim_vec = [sample_struct.primary_sample_flag];

close all
cmap = brewermap(9,'Set2');
samp_fig = figure;
hold on
scatter(region_pca_array(:,1)/1e4,dateNumVec/1e5,'MarkerFaceColor',cmap(8,:),'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2,'MarkerEdgeColor','k')
% scatter(samp_pos_vec,samp_date_vec,'MarkerFaceColor',cmap(3,:),'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2,'MarkerEdgeColor','k')
scatter(samp_pos_vec(samp_prim_vec==1),samp_date_vec(samp_prim_vec==1),'d','MarkerFaceColor',cmap(2,:),'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.2,'MarkerEdgeColor','k')

grid on
xlabel('repo position (au)')
xlabel('repo date (au)')
set(gca,'Fontsize',14)
legend('all repositories','sampled repositories','Location','southeast')

saveas(samp_fig,'../fig/sample_set_dispersion.png')

% [~, si] = sortrows([master_id_vec dateNumVec]);
% x_mean_vec_sorted = x_mean_vec(si);
% y_mean_vec_sorted = y_mean_vec(si);
% [~, dist_score_array] = pca([x_mean_vec_sorted y_mean_vec_sorted]);
% 
% repoCellSorted = repoCell(si);
% regionCellSorted = regionCell(si);
% master_id_vec_sorted = master_id_vec(si);
% subIDVec = subIDVec(si);
% dateNumVec = dateNumVec(si);
% NameIndexCleanSorted = NameIndexClean(si);
% NameIndexRawSorted = NameIndexRaw1(si);
% dateStringCellSorted = dateStringCell(si);
% 
% % Get counts by repo
% [iu_rep, map_to, map_from] = unique(master_id_vec_sorted);
% 
% rep_count_vec_long = NaN(size(map_from));
% rep_count_vec = NaN(size(iu_rep));
% date_count_vec_long = NaN(size(map_from));
% date_diff_vec = NaN(size(iu_rep));
% 
% for i = 1:length(iu_rep)
%   
%     % number of occurrences
%     rep_count_vec_long(map_from==i) = sum(map_from==i);
%     rep_count_vec(i) = sum(map_from==i);
%     
%     if sum(map_from==i) > 1
%         date_diff_vec(i) = median(diff(dateNumVec(map_from==i)));
%         date_count_vec_long(map_from==i) = median(diff(dateNumVec(map_from==i)));
%     end
% end
% 
% %% Get frequency counts for distribution in space and time
% min_reps = 2; % minimum number of times a repo can appear in the dataset
% max_med_date_dist = 180; % maximum median distance between consecutive observations
% nBins = 10;
% 
% % apply filter
% date_rep_filter = find(rep_count_vec >= min_reps & date_diff_vec <= max_med_date_dist);
% 
% % get empircal distribution of counts
% dateNumVecTrunc = dateNumVec(map_to);
% posVecTrunc = dist_score_array(map_to);
% [count_array,~,~,xID,yID] = histcounts2(dateNumVecTrunc,posVecTrunc,nBins);
% 
% % draw random sample weighted by empirical counts
% linIndex = sub2ind(size(count_array),xID,yID);
% rng(421);
% region_index_list = randsample(date_rep_filter,1e4*length(date_rep_filter),true,1./count_array(linIndex(date_rep_filter)));
% 
% % enforce uniqueness
% region_index_list_u = unique(region_index_list,'stable');
% region_id_list_u = iu_rep(region_index_list_u);
% 
% % map back to original dimensions
% % region_name_cell = {}; 
% remapping_vec = [];
% for i = 1:length(region_id_list_u)
%     remapping_vec = [remapping_vec find(map_from==region_id_list_u(i))'];
% %     region_name_cell = [region_name_cell{:} NameIndexCleanSorted(map_from==region_id_list_u(i))'];
% end    
% 
% %% initialize structure
% region_sampling_table = struct;
% 
% master_id_vec_orig = master_id_vec_sorted(remapping_vec);
% [~,~,region_sampling_table.sample_order_vec] = unique(master_id_vec_orig,'stable');
% 
% region_sampling_table.region_name_clean = NameIndexCleanSorted(remapping_vec);
% region_sampling_table.region_name_orig = NameIndexRawSorted(remapping_vec);
% % region_sampling_table.date_num = dateNumVec(remapping_vec);
% region_sampling_table.date_string = dateStringCell(remapping_vec);
% 
% 
% region_sampling_table = struct2table(region_sampling_table);
% 
% % save
% writetable(region_sampling_table,[OutDir 'region_sampling_list.csv']);
% 

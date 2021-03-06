% Script to radnomly select subset of repositories to segment in detail for
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
% [unique_regions, ~, region_id_index] = unique(regionCell);
[region_id_index, ~, master_id_vec] = unique(repoCell);

% Rearrange vectors and identify regions with multiple entries
[~, si] = sortrows([master_id_vec dateNumVec]);
x_mean_vec_sorted = x_mean_vec(si);
y_mean_vec_sorted = y_mean_vec(si);
[~, dist_score_array] = pca([x_mean_vec_sorted y_mean_vec_sorted]);

repoCellSorted = repoCell(si);
regionCellSorted = regionCell(si);
master_id_vec_sorted = master_id_vec(si);
subIDVec = subIDVec(si);
dateNumVec = dateNumVec(si);
NameIndexCleanSorted = NameIndexClean(si);
NameIndexRawSorted = NameIndexRaw1(si);
dateStringCellSorted = dateStringCell(si);

% Get counts by repo
[iu_rep, map_to, map_from] = unique(master_id_vec_sorted);

rep_count_vec_long = NaN(size(map_from));
rep_count_vec = NaN(size(iu_rep));
date_count_vec_long = NaN(size(map_from));
date_diff_vec = NaN(size(iu_rep));

for i = 1:length(iu_rep)
  
    % number of occurrences
    rep_count_vec_long(map_from==i) = sum(map_from==i);
    rep_count_vec(i) = sum(map_from==i);
    
    if sum(map_from==i) > 1
        date_diff_vec(i) = median(diff(dateNumVec(map_from==i)));
        date_count_vec_long(map_from==i) = median(diff(dateNumVec(map_from==i)));
    end
end

%% Get frequency counts for distribution in space and time
min_reps = 2; % minimum number of times a repo can appear in the dataset
max_med_date_dist = 180; % maximum median distance between consecutive observations
nBins = 10;

% apply filter
date_rep_filter = find(rep_count_vec >= min_reps & date_diff_vec <= max_med_date_dist);

% get empircal distribution of counts
dateNumVecTrunc = dateNumVec(map_to);
posVecTrunc = dist_score_array(map_to);
[count_array,~,~,xID,yID] = histcounts2(dateNumVecTrunc,posVecTrunc,nBins);

% draw random sample weighted by empirical counts
linIndex = sub2ind(size(count_array),xID,yID);
rng(421);
region_index_list = randsample(date_rep_filter,1e4*length(date_rep_filter),true,1./count_array(linIndex(date_rep_filter)));

% enforce uniqueness
region_index_list_u = unique(region_index_list,'stable');
region_id_list_u = iu_rep(region_index_list_u);

% map back to original dimensions
% region_name_cell = {}; 
remapping_vec = [];
for i = 1:length(region_id_list_u)
    remapping_vec = [remapping_vec find(map_from==region_id_list_u(i))'];
%     region_name_cell = [region_name_cell{:} NameIndexCleanSorted(map_from==region_id_list_u(i))'];
end    

%% initialize structure
region_sampling_table = struct;

master_id_vec_orig = master_id_vec_sorted(remapping_vec);
[~,~,region_sampling_table.sample_order_vec] = unique(master_id_vec_orig,'stable');

region_sampling_table.region_name_clean = NameIndexCleanSorted(remapping_vec);
region_sampling_table.region_name_orig = NameIndexRawSorted(remapping_vec);
% region_sampling_table.date_num = dateNumVec(remapping_vec);
region_sampling_table.date_string = dateStringCell(remapping_vec);


region_sampling_table = struct2table(region_sampling_table);

% save
writetable(region_sampling_table,[OutDir 'region_sampling_list.csv']);


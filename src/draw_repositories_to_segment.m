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

%% First, clean up location names
refTableClean = refTable;
[NameIndexRaw, ia, ic] = unique(refTable.Name);
NameIndexClean = NameIndexRaw;

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

%% initialize new structure
NameIndexClean = unique(NameIndexClean);
RepoTable.RepoNames = NameIndexClean;

% now parse the names
regionCell = cell(size(NameIndexClean));
dateVec = datetime(zeros(size(NameIndexClean)),0,0);
subIDVec = NaN(size(NameIndexClean));
uniqueIDVec = (1:length(NameIndexClean))';

for i = 1:length(NameIndexClean)
    name = NameIndexClean{i};
    % find underscores   
    underscores = strfind(name,'_');
    % get region name
    regionCell{i} = name(1:underscores(1)-1);
    % extract date
    date_string = name(underscores(1)+1:underscores(2)-1);
    dt = datetime(date_string,'InputFormat','MMddyyyy');
    dateVec(i) = dt;
    % repo sub ID
    subIDVec(i) = str2double(name(underscores(2)+1:underscores(3)-1));    
end    

% Now generate vector of region IDs and a master ID vector
dateNumVec = datenum(dateVec);
[unique_regions, ~, region_id_index] = unique(regionCell);
[region_id_array, ~, master_id_vec] = unique([region_id_index, subIDVec], 'rows');

%% Rearrange vectors and identify regions with multiple entries
[~, si] = sortrows([master_id_vec dateNumVec]);
master_id_vec_sorted = master_id_vec(si);
region_id_index = region_id_index(si);
subIDVec = subIDVec(si);
dateNumVec = dateNumVec(si);

[iu, ~, ix] = unique(master_id_vec_sorted);
count_vec = NaN(size(ix));
for i = 1:length(iu)
    count_vec(ix==i) = sum(ix==i);
end
%%
    

rm_filter = infoTable.plateau_area./infoTable.base_area<0.2; % remove cases when very little plateua area was observed
h_v_a = fitlm(infoTable.base_area(~rm_filter), infoTable.max_height(~rm_filter));
coeff_vec = [h_v_a.Coefficients.Estimate' h_v_a.RMSE];

% extract vectors
x_vec = refTable.POINT_X;
y_vec = refTable.POINT_Y;
a_vec = refTable.area;
o_vec = refTable.ORIG_FID;

% generate unique area vec
region_index = unique(o_vec);
a_vec_u = NaN(size(region_index));
for r = 1:length(region_index)
  a_vec_u(r) = mean(a_vec(o_vec==region_index(r)));
end
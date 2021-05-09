clear
close all

% set paths
RawDataPath = '../dat/';
ProcessedDataPath = '../out/';

% load flux and capacity datasets
load([ProcessedDataPath 'master_capacity_struct.mat'],'master_capacity_struct')
load([ProcessedDataPath 'master_flux_struct.mat'],'master_flux_struct')
load([ProcessedDataPath 'height_fit_struct.mat'],'height_fit_struct')
angle_rep_vec = [height_fit_struct.angle_rep];
%% First let's write basic volume flux and volume capacity results to csv
%array2struct(zeros(1,8),'VariableNames',{'ORIG_FID','repose_angle','vol_est','vol_ste'});
total_vol_array = NaN(length(angle_rep_vec)*length(master_capacity_struct),3);
iter = 1;
for f = 1:length(master_capacity_struct)
    repo_id_index = master_capacity_struct(f).repo_id_index;    
    
    % calculate means and bootstrap standard error
    mean_vol_array = nanmean(master_capacity_struct(f).vol_capacity_array,1);
    ste_vol_array = nanstd(master_capacity_struct(f).vol_capacity_array,[],1);
    
    % generate table
    capacity_array_out = NaN(length(angle_rep_vec)*length(repo_id_index),3);  
    start_i = 1;
    for a = 1:size(mean_vol_array,3)
        last_i = start_i + length(repo_id_index) - 1;
        capacity_array_out(start_i:last_i,1) = repelem(angle_rep_vec(a),length(repo_id_index));
        capacity_array_out(start_i:last_i,2) = mean_vol_array(1,:,a);
        capacity_array_out(start_i:last_i,3) = ste_vol_array(1,:,a);
        
        start_i = last_i + 1;
        
        % get total colume stats
         total_vol_array(iter,1) = length(repo_id_index);
         total_vol_array(iter,2) = angle_rep_vec(a);
         total_vol_array(iter,3) = nanmean(sum(master_capacity_struct(f).vol_capacity_array(:,:,a),2));
         total_vol_array(iter,4) = nanstd(sum(master_capacity_struct(f).vol_capacity_array(:,:,a),2));
         iter = iter + 1;
    end
    % convert to table
    capacity_table = array2table(capacity_array_out,'VariableNames',{'repose_angle','vol_est','vol_ste'});
    capacity_table.ORIG_FID = repmat(repo_id_index,3,1);
    capacity_table = [capacity_table(:,end) capacity_table(:,1:end-1)];
    
    % write to file
    source_name = master_capacity_struct(f).source_file;
%     source_name = source_name(1:end-4);
    writetable(capacity_table,['ind_repo_capacity_' source_name]); 
end

total_vol_table = array2table(total_vol_array,'VariableNames',{'number_of_repos','repose_angle','total_vol_est','total_vol_ste'});
total_vol_table.source_file = repelem({master_capacity_struct.source_file},length(angle_rep_vec))';
total_vol_table = [total_vol_table(:,end) total_vol_table(:,1:end-1)];
writetable(total_vol_table,'total_capacity_estimates.csv'); 
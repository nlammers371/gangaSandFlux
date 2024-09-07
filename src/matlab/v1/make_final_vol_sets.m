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
    writetable(capacity_table,[ProcessedDataPath 'ind_repo_capacity_' source_name]); 
end

total_vol_table = array2table(total_vol_array,'VariableNames',{'number_of_repos','repose_angle','total_vol_est','total_vol_ste'});
total_vol_table.source_file = repelem({master_capacity_struct.source_file},length(angle_rep_vec))';
total_vol_table = [total_vol_table(:,end) total_vol_table(:,1:end-1)];
writetable(total_vol_table,[ProcessedDataPath 'total_capacity_estimates.csv']); 

%% Now do same for for volume flux estimates
repo_vol_array = [];
date_index = master_flux_struct(1).date_index;
repo_id_index = master_flux_struct(1).repo_id_index;
% generate strings to use for repo variable names
repo_name_cell = {};
for r = 1:1:length(repo_id_index)
    repo_name_cell(r) = {['repo_' sprintf('%02d',r) '_vol_mean']};
    repo_name_cell(length(repo_id_index)+r) = {['repo_' sprintf('%02d',r) '_vol_ste']};
end    
% generate date strings
date_name_cell = {};
for d = 1:length(date_index)
    date_name_cell(d) = {datestr(date_index(d))};
end
dt_val = date_index(end)-date_index(1);
dt_vec = diff(date_index);
[min_dt ,min_dt_i] = min(dt_vec);

volum_flux_array = NaN(length(angle_rep_vec),12);
for a = 1:length(master_flux_struct)
    
    vol_array = master_flux_struct(a).vol_flux_array;
    repo_capacity_array = permute(master_flux_struct(a).vol_fp_array,[3 2 1]);
    vol_array_rel = vol_array./repo_capacity_array * 100;
    repo_capacity_mean = nanmean(repo_capacity_array,3)';
    % calculate and save stats for individual repos
    vol_mean_array = reshape(nanmean(vol_array,3),[],1);
    vol_ste_array = reshape(nanstd(vol_array,[],3),[],1);
    vol_array = [vol_mean_array vol_ste_array repelem(repo_capacity_mean,length(date_index))];
    vol_table = array2table(vol_array,'VariableNames', {'current_vol_est', 'current_vol_ste','total_capacity'});
    dates_long = repmat(date_index,length(repo_id_index),1);
    vol_table.date_num = dates_long;
    vol_table.date_string = repmat(date_name_cell',length(repo_id_index),1);
    vol_table.repo_id = repelem(repo_id_index,length(date_index))';
    vol_table = [vol_table(:,end) vol_table(:,end-2:end-1) vol_table(:,1:end-3)];
    
    writetable(vol_table,[ProcessedDataPath 'repo_volume_vs_time_angle_' num2str(angle_rep_vec(a)) '.csv']); 
    
    % calculate population stats
    
    % flux in   
    % calculate cumulative flux metrics
    vol_flux_array_rel = diff(vol_array_rel,1,1);    
    vol_flux_in_array_rel = vol_flux_array_rel;
    vol_flux_in_array_rel(vol_flux_in_array_rel<0) = 0;
    
    % get stats for total observed cohort
    vol_flux_in_total_mean = nanmean(nanmean(vol_flux_in_array_rel,2),3);
    vol_flux_in_total_ste = nanstd(nanmean(vol_flux_in_array_rel,2),[],3);
    volum_flux_array(a,1) = sum(vol_flux_in_total_mean);
    volum_flux_array(a,2) = sum(vol_flux_in_total_ste);
    
    % flux out
    vol_flux_out_array_rel = vol_flux_array_rel;
    vol_flux_out_array_rel(vol_flux_out_array_rel>0) = 0;
    
    % get stats for total observed cohort
    vol_flux_out_total_mean = nanmean(nanmean(vol_flux_out_array_rel,2),3);
    vol_flux_out_total_ste = nanstd(nanmean(vol_flux_out_array_rel,2),[],3);
    volum_flux_array(a,3) = sum(vol_flux_out_total_mean);
    volum_flux_array(a,4) = sum(vol_flux_out_total_ste);
    
    % replicate these numbers, but divided by time to give a rate
    volum_flux_array(a,5:8) = volum_flux_array(a,1:4)/dt_val;
    
    % now calculate a less conservative flux estimate    
    d_out_dt = nanmean(nanmean(vol_flux_out_array_rel./dt_vec,2),3);
    d_out_dt_ste = nanstd(nanmean(vol_flux_out_array_rel./dt_vec,2),[],3);
    d_in_dt = nanmean(nanmean(vol_flux_in_array_rel./dt_vec,2),3);
    d_in_dt_ste = nanstd(nanmean(vol_flux_in_array_rel./dt_vec,2),[],3);
        
    volum_flux_array(a,9:10) = [d_in_dt(min_dt_i) d_in_dt_ste(min_dt_i)];
    volum_flux_array(a,11:12) = [d_out_dt(min_dt_i) d_out_dt_ste(min_dt_i)];
      
end    

% write total flux numbers to file
flux_table = array2table(volum_flux_array,'VariableNames', {'observed_flux_tot_in', 'observed_flux_tot_in_ste','observed_flux_tot_out', 'observed_flux_tot_out_ste',...
                                                                  'observed_flux_rate_in', 'observed_flux_rate_in_ste','observed_flux_rate_out', 'observed_flux_rate_out_ste',...
                                                                  'est_flux_rate_in', 'est_flux_rate_in_ste','est_flux_rate_out', 'est_flux_rate_out_ste'});
flux_table.repose_angle = angle_rep_vec';                                                                
flux_table = [flux_table(:,end) flux_table(:,1:end-1)];
writetable(flux_table,[ProcessedDataPath 'repo_flux_estimates.csv']); 
                                                                
%% Put the pieces together to estimate total sand flux
total_flux_table = total_vol_table;
a_index_vol = total_flux_table.repose_angle;
a_index_flux = flux_table.repose_angle;

% initialize fields
% out
total_flux_table.flux_out_obs = NaN(size(a_index_vol));
total_flux_table.flux_out_obs_ste = NaN(size(a_index_vol));
total_flux_table.flux_out_est = NaN(size(a_index_vol));
total_flux_table.flux_out_est_ste = NaN(size(a_index_vol));
% in
total_flux_table.flux_in_obs = NaN(size(a_index_vol));
total_flux_table.flux_in_obs_ste = NaN(size(a_index_vol));
total_flux_table.flux_in_est = NaN(size(a_index_vol));
total_flux_table.flux_in_est_ste = NaN(size(a_index_vol));

% layer on the flux info 
for a = 1:length(angle_rep_vec)
    a_filter_vol = a_index_vol==angle_rep_vec(a);
    a_filter_flux = a_index_flux==angle_rep_vec(a);
    
    % flux out
    % calculate expected volume per day
    total_flux_table.flux_out_obs(a_filter_vol) = total_flux_table.total_vol_est(a_filter_vol)*flux_table.observed_flux_rate_out(a_filter_flux)/100;
    total_flux_table.flux_out_est(a_filter_vol) = total_flux_table.total_vol_est(a_filter_vol)*flux_table.est_flux_rate_out(a_filter_flux)/100;
    % use linear error propagation to estimate error
    total_flux_table.flux_out_obs_ste(a_filter_vol) = sqrt((total_flux_table.total_vol_est(a_filter_vol)*flux_table.observed_flux_rate_out_ste(a_filter_flux)/100).^2 + ...
                                                           (total_flux_table.total_vol_ste(a_filter_vol)*flux_table.observed_flux_rate_out(a_filter_flux)/100).^2);
                                                         
    total_flux_table.flux_out_est_ste(a_filter_vol) = sqrt((total_flux_table.total_vol_est(a_filter_vol)*flux_table.est_flux_rate_out_ste(a_filter_flux)/100).^2 + ...
                                                           (total_flux_table.total_vol_ste(a_filter_vol)*flux_table.est_flux_rate_out(a_filter_flux)/100).^2);
    
    % flux in
    total_flux_table.flux_in_obs(a_filter_vol) = total_flux_table.total_vol_est(a_filter_vol)*flux_table.observed_flux_rate_in(a_filter_flux)/100;
    total_flux_table.flux_in_est(a_filter_vol) = total_flux_table.total_vol_est(a_filter_vol)*flux_table.est_flux_rate_in(a_filter_flux)/100;
    
    % use linear error propagation to estimate error
    total_flux_table.flux_in_obs_ste(a_filter_vol) = sqrt((total_flux_table.total_vol_est(a_filter_vol)*flux_table.observed_flux_rate_in_ste(a_filter_flux)/100).^2 + ...
                                                           (total_flux_table.total_vol_ste(a_filter_vol)*flux_table.observed_flux_rate_in(a_filter_flux)/100).^2);
                                                         
    total_flux_table.flux_in_est_ste(a_filter_vol) = sqrt((total_flux_table.total_vol_est(a_filter_vol)*flux_table.est_flux_rate_in_ste(a_filter_flux)/100).^2 + ...
                                                           (total_flux_table.total_vol_ste(a_filter_vol)*flux_table.est_flux_rate_in(a_filter_flux)/100).^2);
end    

writetable(total_flux_table,[ProcessedDataPath 'total_daily_flux_estimates.csv']); 


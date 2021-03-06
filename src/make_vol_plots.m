% Script to make basic plots of SA, volume, and flux over time
clear
close all

DataPath = '../dat/20220214/';
FileName = 'sbapts_02142022putm.csv';
OutDir = '../out/20220214/';

FigDir = '../fig/20220214/';
mkdir(FigDir)

% load table
refTable = readtable([DataPath FileName]);
volTable = readtable([OutDir 'repo_volume_estimates.csv']);

problem_repo_names = {'n_Bheramara_01172011_560_b','n_Bheramara_05142011_369_b','n_Bheramara_11042011_243_b'};

% Generate summary stat table
nBoots = 100;
region_date_array = unique([volTable.region_id volTable.date_num],'rows');

name_flags = ismember(volTable.names_raw,problem_repo_names);

total_vol_mean = NaN(size(region_date_array,1),1);
total_vol_ste = NaN(size(region_date_array,1),1);
mean_vol_mean = NaN(size(region_date_array,1),1);
mean_vol_ste = NaN(size(region_date_array,1),1);
total_sa_array = NaN(size(region_date_array,1),1);
mean_sa_array = NaN(size(region_date_array,1),1);
n_repo_array = NaN(size(region_date_array,1),1);
mean_x_array = NaN(size(region_date_array,1),1);
mean_y_array = NaN(size(region_date_array,1),1);
date_vec = NaN(size(region_date_array,1),1);
region_vec = NaN(size(region_date_array,1),1);

for r = 1:size(region_date_array,1)
  
    % get filter
    r_filter = volTable.region_id==region_date_array(r,1) & volTable.date_num==region_date_array(r,2) & ~name_flags;
    
    % extract vol data
    vol_mean_vec = volTable.vol_capacity_mean(r_filter);
    vol_ste_vec = volTable.vol_capacity_ste(r_filter);
    
    % generate random samples
    vol_mean_array = normrnd(repmat(vol_mean_vec,1,nBoots),repmat(vol_ste_vec,1,nBoots));
    
    % save volume results
    total_vol_mean(r) = mean(sum(vol_mean_array,1));
    total_vol_ste(r) = std(sum(vol_mean_array,1));
    mean_vol_mean(r) = mean(mean(vol_mean_array,1));
    total_vol_ste(r) = std(mean(vol_mean_array,1));
    
    % save other results
    total_sa_array(r) = sum(volTable.sa_capacity(r_filter));
    mean_sa_array(r) = mean(volTable.sa_capacity(r_filter));
    n_repo_array(r) = sum(r_filter);
    date_vec(r) = region_date_array(r,2);
    region_vec(r) = region_date_array(r,1);
    mean_x_array(r) = mean(volTable.mean_x_pos(r_filter));
    mean_y_array(r) = mean(volTable.mean_y_pos(r_filter));
    
end    

% Generate region-specific volume trend plots
close all
region_index = unique(region_vec);
region_trend_dir = [FigDir 'region_trends' filesep]; 
mkdir(region_trend_dir);

p_param_range = linspace(1e-5,1e-7);

date_axis = min(date_vec):max(date_vec);
vol_mean_array = NaN(length(date_axis),length(region_index));
vol_ste_array = NaN(length(date_axis),length(region_index));

for r = 1:length(region_index)
       
    
    region_fig = figure('Visible','off');  
    hold on
    cmap = brewermap(8,'Set2');
    
    % extract info
    dt = date_vec(region_vec==region_index(r));    
    tv = total_vol_mean(region_vec==region_index(r));    
    tve = total_vol_ste(region_vec==region_index(r));
            
    % calculate smoothing spline
    dt_indices = find(dt(1)==date_axis):find(dt(end)==date_axis);
    
    temp_array = NaN(length(dt_indices),nBoots);
    
    if length(dt_indices)>1        
        for n = 1:nBoots
            tv_temp = normrnd(tv,tve);
            temp_array(:,n) = csaps(dt,tv_temp,randsample(p_param_range,1,true),date_axis(dt_indices));        
        end
    else
        temp_array(:,n) = tv;
    end
    dt_plot = datetime(date_axis(dt_indices),'ConvertFrom','datenum');
    ub95 = prctile(temp_array,95,2);
    lb05 = prctile(temp_array,5,2);
    tv_trend = nanmean(temp_array,2);
    
    fill([dt_plot fliplr(dt_plot)], [ub95' fliplr(lb05')],'k','FaceAlpha',0.3,'EdgeAlpha',0)
    plot(dt_plot,tv_trend,'Color','k','LineWidth',1.5)
    scatter(datetime(dt,'ConvertFrom','datenum'),tv,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k')
 
    grid on
    set(gca,'FontSize',14)
    
    ylabel('total capacity in region (m^3)')
    xlabel('date')    
    saveas(region_fig,[region_trend_dir 'capacity_trend_region' sprintf('%03d',region_index(r)) '.png'])

    % record
    vol_mean_array(dt_indices,r) = tv_trend;
    vol_ste_array(dt_indices,r) = nanstd(temp_array,[],2);
end    

%% Compile volume and flux statistics

close all
vol_mean_array(vol_mean_array<0) = 0;
vol_mean_array_filled = vol_mean_array;
vol_ste_array_filled = vol_ste_array;
for i = 1:size(vol_mean_array_filled,2)
    start_i = find(~isnan(vol_mean_array(:,i)),1);
    last_i = find(~isnan(vol_mean_array(:,i)),1,'last');
    
    vol_mean_array_filled(1:start_i-1,i) = 0;
    vol_mean_array_filled(last_i+1:end,i) = vol_mean_array_filled(last_i,i);
    
    vol_ste_array_filled(1:start_i-1,i) = 0;
    vol_ste_array_filled(last_i+1:end,i) = vol_ste_array_filled(last_i,i);
end
% vol_mean_array_filled(vol_mean_array_filled<0) = 0;

vol_total_array_boot = NaN(size(vol_mean_array_filled,1),nBoots);
for n = 1:nBoots
    vol_total_array_boot(:,n) = sum(normrnd(vol_mean_array_filled,vol_ste_array_filled),2);
end

%% Make figures
close all
date_axis_plot = datetime(date_axis,'ConvertFrom','datenum');

% plot all repo volumes together
region_fig = figure;
hold on
for i = 1:size(vol_mean_array,2)
    plot(date_axis_plot, vol_mean_array(:,i));
end

xlabel('date')
grid on
set(gca,'FontSize',14)    
ylabel('total capacity by region (m^3)')
xlabel('date')

xtickangle(-45)
region_fig.Renderer='Painters';

saveas(region_fig,[FigDir 'region_capacity_plot.png'])
saveas(region_fig,[FigDir 'region_capacity_plot.pdf'])

% plot total volume across all regions
total_vol_mean = nanmean(vol_total_array_boot,2);
total_vol_ste = nanstd(vol_total_array_boot,[],2);
ub = total_vol_mean+total_vol_ste;
lb = total_vol_mean-total_vol_ste;

%%
total_fig = figure;
area(date_axis_plot,imgaussfilt(total_vol_mean,5),'FaceColor',cmap(7,:),'LineWidth',1,'EdgeColor','k')
% fill([date_axis_plot fliplr(date_axis_plot)], [ub' fliplr(lb')],'k','FaceAlpha',0.3,'EdgeAlpha',0)
% plot(date_axis_plot,total_vol_mean,'Color','k','LineWidth',2)

xlabel('date')
grid on
set(gca,'FontSize',14)    
ylabel('total capacity (m^3)')
xlabel('date')

saveas(total_fig,[FigDir 'total_capacity_plot.png'])
saveas(total_fig,[FigDir 'total_capacity_plot.pdf'])

% plot approximate daily flux
flux_pct_vec = [0.01 0.02 0.035];
c_vec = [2 3 5];
for f = 1:length(flux_pct_vec)
  
    flux_fig = figure;    
    area(date_axis_plot,imgaussfilt(total_vol_mean*flux_pct_vec(f),5),'FaceColor',cmap(c_vec(f),:),'LineWidth',1,'EdgeColor','k')

    xlabel('date')
    grid on
    set(gca,'FontSize',14)    
    ylabel('total daily flux out (m^3)')
    xlabel('date')
    yl = get(gca,'YLim');
    ylim([0 yl(2)]);
    saveas(flux_fig,[FigDir 'total_flux_plot_f' sprintf('%03d',round(100*flux_pct_vec(f),0)) '.png'])
    saveas(flux_fig,[FigDir 'total_flux_plot_f' sprintf('%03d',round(100*flux_pct_vec(f),0)) '.pdf'])
end    

%% attempt to plot cumulative 
month_vec = month(date_axis_plot);
plot_months = [11:12 1:4];
for f = 1:length(flux_pct_vec)
    flux_vec = total_vol_mean*flux_pct_vec(f);
    flux_vec(~ismember(month_vec,plot_months)) = zeros;
    cum_flux_vec = cumsum(flux_vec);
    flux_reset = interp1(date_axis(~ismember(month_vec,plot_months))',cum_flux_vec(~ismember(month_vec,plot_months)),date_axis','previous');
   
    clum_flux_reset = cum_flux_vec-flux_reset;
    
    flux_fig = figure;    
    area(date_axis_plot,clum_flux_reset,'FaceColor',cmap(c_vec(f),:),'LineWidth',1,'EdgeColor','k')

    xlabel('date')
    grid on
    set(gca,'FontSize',14)    
    ylabel('cumulative flux out (m^3)')
    xlabel('date')
    yl = get(gca,'YLim');
    ylim([0 yl(2)]);
    xtickangle(-45)
    saveas(flux_fig,[FigDir 'cycle_flux_plot_f' sprintf('%03d',round(100*flux_pct_vec(f),0)) '.png'])
    saveas(flux_fig,[FigDir 'cycle_flux_plot_f' sprintf('%03d',round(100*flux_pct_vec(f),0)) '.pdf'])
end  



% cmap = flipud(brewermap(8,'Spectral'));
% 
% vol_plot = volTable.vol_capacity_mean; 
% vol_99 = prctile(vol_plot,99.5);
% vol_plot(vol_plot>vol_99) = vol_99;
% vol_plot = vol_plot / vol_99 * 200;
% 
% scatter_fig = figure;
% colormap(cmap)
% scatter(volTable.mean_x_pos,volTable.mean_y_pos,vol_plot, volTable.date_num...
%         ,'filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
%       
% grid on
% xlabel('x position')
% ylabel('y position')
% saveas(scatter_fig,[FigDir 'repo_scatter.png'])
clear
close all

ReadPath = '../out/';
FigDir = '../fig/';
mkdir(FigDir)

load([ReadPath 'master_flux_struct.mat'],'master_flux_struct')
load([ReadPath 'height_fit_struct.mat'],'height_fit_struct')

%% plot cumulative observed flux in and out
angle_rep_index = 2;

% load variables into workspace
fnames = fieldnames(master_flux_struct);
for f = 1:length(fnames)
    eval([fnames{f} ' = master_flux_struct(angle_rep_index).(fnames{f});'])
end

% out
vol_flux_out_array = cat(1,zeros(1,length(repo_id_index),n_boots),diff(vol_flux_array,1,1));
vol_flux_out_array(vol_flux_out_array>0) = 0;
vol_flux_out_array_mean = nanmean(vol_flux_out_array,3)./nanmean(vol_fp_array)*100;
vol_flux_out_array_ste = nanstd(vol_flux_out_array,[],3)./nanmean(vol_fp_array)*100;
cum_vol_flux_out_array = cumsum(vol_flux_out_array_mean);
cum_vol_flux_out_array_ste = cumsum(vol_flux_out_array_ste);

% calculate cumulative flux metrics
vol_flux_out_array_rel = vol_flux_out_array./nanmean(vol_fp_array)*100;
% get stats for total observed cohort
vol_flux_out_total_mean = nanmean(nanmean(vol_flux_out_array_rel,2),3);
vol_flux_out_total_ste = nanstd(nanmean(vol_flux_out_array_rel,2),[],3);
cum_vol_flux_out_total = cumsum(vol_flux_out_total_mean);
cum_vol_flux_out_total_ste = cumsum(vol_flux_out_total_ste);

% in
vol_flux_in_array = cat(1,zeros(1,length(repo_id_index),n_boots),diff(vol_flux_array,1,1));
vol_flux_in_array(vol_flux_in_array<0) = 0;

vol_flux_in_array_rel = vol_flux_in_array./nanmean(vol_fp_array)*100;

% get stats for individual replicates
vol_flux_in_array_mean = nanmean(vol_flux_in_array,3)./nanmean(vol_fp_array)*100;
vol_flux_in_array_ste = nanstd(vol_flux_in_array,[],3)./nanmean(vol_fp_array)*100;
cum_vol_flux_in_array = cumsum(vol_flux_in_array_mean);
cum_vol_flux_in_array_ste = cumsum(vol_flux_in_array_ste);

% get stats for total observed cohort
vol_flux_in_total_mean = nanmean(nanmean(vol_flux_in_array_rel,2),3);
vol_flux_in_total_ste = nanstd(nanmean(vol_flux_in_array_rel,2),[],3);
cum_vol_flux_in_total = cumsum(vol_flux_in_total_mean);
cum_vol_flux_in_total_ste = cumsum(vol_flux_in_total_ste);

close all
date_plot = date_index'-date_index(1);

flux_fig = figure;
cmap_in = brewermap(length(repo_id_index)+1,'Reds'); 
cmap_out = brewermap(length(repo_id_index)+1,'Blues'); 

hold on
for r = 1:length(repo_id_index)    
    ub_vec = cum_vol_flux_out_array(:,r) + cum_vol_flux_out_array_ste(:,r);
    lb_vec = cum_vol_flux_out_array(:,r) - cum_vol_flux_out_array_ste(:,r);
    
    fill([date_plot fliplr(date_plot)],[ub_vec',fliplr(lb_vec')],cmap_out(r,:),'FaceAlpha',0.2,'EdgeAlpha',0.3)
    plot(date_index-date_index(1),cum_vol_flux_out_array(:,r),'Color',cmap_out(r,:),'LineWidth',1.5)
    scatter(date_index-date_index(1),cum_vol_flux_out_array(:,r),'MarkerFaceColor',cmap_out(r,:),'MarkerEdgeColor','k')    
    
    
    ub_vec = cum_vol_flux_in_array(:,r) + cum_vol_flux_in_array_ste(:,r);
    lb_vec = cum_vol_flux_in_array(:,r) - cum_vol_flux_in_array_ste(:,r);
    
    fill([date_plot fliplr(date_plot)],[ub_vec',fliplr(lb_vec')],cmap_in(r,:),'FaceAlpha',0.2,'EdgeAlpha',0.3)
    plot(date_index-date_index(1),cum_vol_flux_in_array(:,r),'Color',cmap_in(r,:),'LineWidth',1.5)  
    scatter(date_index-date_index(1),cum_vol_flux_in_array(:,r),'MarkerFaceColor',cmap_in(r,:),'MarkerEdgeColor','k')
%     scatter(date_index,vol_flux_array(:,r)./vol_flux_array(1,r),'MarkerFaceColor',cmap(r,:),'MarkerEdgeColor','k')
end

grid on
ylim([-300 300])
set(gca,'Fontsize',12)
xlabel('days since first observation')
ylabel('cumulative change (% repo capacity)')
saveas(flux_fig,[FigDir 'cumulative_relative_flux.png'])
saveas(flux_fig,[FigDir 'cumulative_relative_flux.pdf'])
  
% same plot with means overlaid
flux_tot_fig = figure;

hold on
for r = 1:length(repo_id_index)    
            
    plot(date_index-date_index(1),cum_vol_flux_out_array(:,r),'Color',[cmap_out(r,:) .2],'LineWidth',1)
               
    plot(date_index-date_index(1),cum_vol_flux_in_array(:,r),'Color',[cmap_in(r,:) .2],'LineWidth',1)      
%     scatter(date_index,vol_flux_array(:,r)./vol_flux_array(1,r),'MarkerFaceColor',cmap(r,:),'MarkerEdgeColor','k')
end

% plot totals
c_ind = floor(length(repo_id_index)/2);
ub_vec_out = cum_vol_flux_out_total + cum_vol_flux_out_total_ste;
lb_vec_out = cum_vol_flux_out_total - cum_vol_flux_out_total_ste;
fill([date_plot fliplr(date_plot)],[ub_vec_out',fliplr(lb_vec_out')],cmap_out(c_ind,:),'FaceAlpha',0.5,'EdgeColor','k','EdgeAlpha',0.3)
plot(date_index-date_index(1),cum_vol_flux_out_total,'Color',[cmap_out(c_ind,:) 1],'LineWidth',2)        
scatter(date_index-date_index(1),cum_vol_flux_out_total,'MarkerFaceColor',cmap_out(c_ind,:),'MarkerEdgeColor','k')

ub_vec_in = cum_vol_flux_in_total + cum_vol_flux_in_total_ste;
lb_vec_in = cum_vol_flux_in_total - cum_vol_flux_in_total_ste;
fill([date_plot fliplr(date_plot)],[ub_vec_in',fliplr(lb_vec_in')],cmap_in(c_ind,:),'FaceAlpha',0.5,'EdgeColor','k','EdgeAlpha',0.3)
plot(date_index-date_index(1),cum_vol_flux_in_total,'Color',[cmap_in(c_ind,:) 1],'LineWidth',2)   
scatter(date_index-date_index(1),cum_vol_flux_in_total,'MarkerFaceColor',cmap_in(c_ind,:),'MarkerEdgeColor','k')

grid on
ylim([-210 210])
set(gca,'Fontsize',12)
xlabel('days since first observation')
ylabel('cumulative change (% repo capacity)')
saveas(flux_tot_fig,[FigDir 'mean_cumulative_relative_flux.png'])
saveas(flux_tot_fig,[FigDir 'mean_cumulative_relative_flux.pdf'])

%% Plot inferred flux rate vs point-over-point separation
dt_vec = diff(date_index);
d_out_dt = nanmean(diff(cum_vol_flux_out_array,1,1)./dt_vec,2);
d_in_dt = nanmean(diff(cum_vol_flux_in_array,1,1)./dt_vec,2);

close all
app_flux_fig = figure;
hold on
scatter(dt_vec,d_in_dt,'MarkerFaceColor',cmap_in(c_ind,:),'MarkerEdgeColor','k')
scatter(dt_vec,d_out_dt,'MarkerFaceColor',cmap_out(c_ind,:),'MarkerEdgeColor','k')

set(gca,'Fontsize',12)
xlabel('time difference between observations')
ylabel('apparent volume flux (% capacity per day)')
grid on

saveas(app_flux_fig,[FigDir 'apparent_flux_vs_dt.png'])
saveas(app_flux_fig,[FigDir 'apparent_flux_vs_dt.pdf'])

max_in_rate = nanmax(d_in_dt);
in_flux_projected = date_plot.*max_in_rate;

max_out_rate = nanmin(d_out_dt);
out_flux_projected = date_plot.*max_out_rate;

projected_flux_fig = figure;
hold on

% out
c_ind = floor(length(repo_id_index)/2);
ub_vec_out = cum_vol_flux_out_total + cum_vol_flux_out_total_ste;
lb_vec_out = cum_vol_flux_out_total - cum_vol_flux_out_total_ste;
fill([date_plot fliplr(date_plot)],[ub_vec_out',fliplr(lb_vec_out')],cmap_out(c_ind,:),'FaceAlpha',0.5,'EdgeColor','k','EdgeAlpha',0.3)
plot(date_index-date_index(1),cum_vol_flux_out_total,'Color',[cmap_out(c_ind,:) 1],'LineWidth',2)        
scatter(date_index-date_index(1),cum_vol_flux_out_total,'MarkerFaceColor',cmap_out(c_ind,:),'MarkerEdgeColor','k')

plot(date_plot,out_flux_projected,'--','Color',[cmap_out(c_ind,:) .6],'LineWidth',2)  

% in
ub_vec_in = cum_vol_flux_in_total + cum_vol_flux_in_total_ste;
lb_vec_in = cum_vol_flux_in_total - cum_vol_flux_in_total_ste;
fill([date_plot fliplr(date_plot)],[ub_vec_in',fliplr(lb_vec_in')],cmap_in(c_ind,:),'FaceAlpha',0.5,'EdgeColor','k','EdgeAlpha',0.3)
plot(date_index-date_index(1),cum_vol_flux_in_total,'Color',[cmap_in(c_ind,:) 1],'LineWidth',2)   
scatter(date_index-date_index(1),cum_vol_flux_in_total,'MarkerFaceColor',cmap_in(c_ind,:),'MarkerEdgeColor','k')

plot(date_plot,in_flux_projected,'--','Color',[cmap_in(c_ind,:) .6],'LineWidth',2)  

grid on
ylim([-450 450])
set(gca,'Fontsize',12)
xlabel('days since first observation')
ylabel('cumulative change (% repo capacity)')
saveas(projected_flux_fig,[FigDir 'mean_cumulative_flux_projected.png'])
saveas(projected_flux_fig,[FigDir 'mean_cumulative_flux_projected.pdf'])


%% %%%%%%%%%%%%%% make iullustrative footprint plots %%%%%%%%%%%%%%

footprint_path = [FigDir 'footprint_maps' filesep];
mkdir(footprint_path)

for r = 1:length(repo_id_index)
    stack = repo_flux_struct(r).height_stack(:,:,1);
    mh = nanmax(stack(:));
    for d = 1:length(date_index)
        
        surf_fig = figure('Visible','off');
        cmap = colormap('pink');
        cmap = cmap(75:end,:);
        colormap(cmap);
        surf(repo_flux_struct(r).height_stack(:,:,d))
        
        xlabel('x distance (meters)')
        ylabel('y distance (meters)')
        zlabel('inferred heigh (meters)')
        % h = colorbar;
        % ylabel(h,'plateau height (meters)')

        zlim([0 1.2*mh])

        set(gca,'Fontsize',10)

        saveas(surf_fig,[FigDir 'pile_reconstruction_region_' num2str(r) '_date_' num2str(d) '.png'])
        
        hm_fig = figure('Visible','off');
        cmap1 = flipud(brewermap([],'Spectral'));    
        colormap(cmap1);
        imagesc(repo_flux_struct(r).height_stack(:,:,d))
        
        xlabel('x distance (meters)')
        ylabel('y distance (meters)')
%         zlabel('inferred heigh (meters)')
        h = colorbar;
        ylabel(h,'pile height (meters)')

        zlim([0 10])

        set(gca,'Fontsize',10)

        saveas(hm_fig,[footprint_path 'pile_heatmap_region_' num2str(r) '_date_' num2str(d) '.png'])
    end
end
clear
close all

ProcessedDataPath = '../out/';
DataPath = '../dat/';
FigDir = '../fig/';
mkdir(FigDir)
RefName = 'Rep_Stats_all.xlsx';

load([ProcessedDataPath 'master_flux_struct.mat'],'master_flux_struct')
load([ProcessedDataPath 'height_fit_struct.mat'],'height_fit_struct')
load([ProcessedDataPath 'master_capacity_struct.mat'],'master_capacity_struct')
infoTable = readtable([DataPath RefName]);
angle_rep_index = 2;
%% Make repo height fit plot 

area_ref_vec = infoTable.base_area;
area_index = linspace(min(area_ref_vec),2*max(area_ref_vec));

close all
height_fig = figure;
cmap = brewermap([],'Set2');
hold on
p1 = plot(area_index,height_fit_struct(angle_rep_index).max_height_vec,'Color','k','LineWidth',1.5);
s1 = scatter(area_ref_vec,infoTable.height_30,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k');

grid on
set(gca,'Fontsize',14)
xlabel('repository surface area (m^2)')
ylabel('estimate height (m)')
legend([s1 p1],'indiividual repos','model fit','Location','southeast')
xlim([0 3e4])
saveas(height_fig,[FigDir 'repo_max_height_fit.png'])
saveas(height_fig,[FigDir 'repo_max_height_fit.pdf'])

%% Make histograms of repo volume for each dataset
bins = linspace(0,6e4,101);
vol_array = NaN(length(master_capacity_struct),length(bins)-1);
% for i = 1:length(master_capacity_struct)
%   hist_fig = figure;
%   cmap = brewermap([],'Set2');
%   vol_vec = nanmean(master_capacity_struct(i).vol_capacity_array(:,:,angle_rep_index));
%   h = histogram(vol_vec,bins,'Normalization','probability','FaceColor',cmap(i,:));
%   vol_array(i,:) = h.Values;
%   xlabel('repository capacity (m^3)')
%   ylabel('probability')
%   set(gca,'Fontsize',14)
%   grid on
%   suffix = master_capacity_struct(i).source_file(1:end-4);
%   saveas(hist_fig,[FigDir 'repo_max_height_fit_' suffix '.png'])
%   saveas(hist_fig,[FigDir 'repo_max_height_fit_' suffix '.pdf'])
% end

% plot together as cumulative distribtutions
% make logarithmic bins
bins_log = logspace(1,5,151);
x_axis = bins_log(1:end-1) + diff(bins_log)/2;

cum_fig = figure;
cmap = brewermap([],'Set2');
hold on
for i = 1:length(master_capacity_struct)
  vol_vec = nanmean(master_capacity_struct(i).vol_capacity_array(:,:,angle_rep_index));
  counts = histcounts(vol_vec,bins_log);
  counts = counts/sum(counts);
  plot(x_axis,cumsum(counts),'Color',cmap(i,:),'LineWidth',2)
end
set(gca,'xScale','log');
set(gca,'Fontsize',14)
grid on
xlabel('repository capacity (m^3)')
ylabel('cumulative probability (eCDF)')
legend('2020','2015','2009-2010','Location','southeast')

saveas(cum_fig,[FigDir 'volume_cdf.png'])
saveas(cum_fig,[FigDir 'volume_cdf.pdf'])

%% make illustrative footprint plots

% load variables into workspace
fnames = fieldnames(master_flux_struct);
for f = 1:length(fnames)
    eval([fnames{f} ' = master_flux_struct(angle_rep_index).(fnames{f});'])
end

footprint_path = [FigDir 'footprint_maps_main' filesep];
mkdir(footprint_path)

plot_ids = [14 15];
close all
for r = plot_ids
    stack = repo_flux_struct(r).height_stack(:,:,1);
    mh = nanmax(stack(:));
    for d = 1:length(date_index)
        
        surf_fig = figure;%('Visible','off');
        cmap = colormap('pink');
        cmap = cmap(75:end,:);
        colormap(cmap);
        height_slice = imgaussfilt(repo_flux_struct(r).height_stack(:,:,d),1);
        surf(height_slice)        
        xlabel('x distance (meters)')
        ylabel('y distance (meters)')
        zlabel('inferred heigh (meters)')
        % h = colorbar;
        % ylabel(h,'plateau height (meters)')

        zlim([0 6])
        xlim([0 size(height_slice,2)])
        ylim([0 size(height_slice,1)])
        view(225,35)
        set(gca,'Fontsize',10)

        saveas(surf_fig,[footprint_path 'pile_reconstruction_region_' num2str(r) '_date_' num2str(d) '.png'])
        
        hm_fig = figure;%('Visible','off');
        cmap1 = flipud(brewermap([],'Spectral'));    
        colormap(cmap);
        imagesc(flipud(repo_flux_struct(r).height_stack(:,:,d)')==0)
        
        xlabel('x distance (meters)')
        ylabel('y distance (meters)')
%         zlabel('inferred heigh (meters)')
%         h = colorbar;
%         ylabel(h,'pile height (meters)')

        zlim([0 10])

        set(gca,'Fontsize',12)

        saveas(hm_fig,[footprint_path 'pile_mask_' num2str(r) '_date_' num2str(d) '.png'])
    end
end

%% Plot an overlay for 15 to show process
fnames = fieldnames(master_flux_struct);
for f = 1:length(fnames)
    eval([fnames{f} ' = master_flux_struct(angle_rep_index).(fnames{f});'])
end

plot_id =  15;
date_ids = [7 1];

close all

stack = repo_flux_struct(plot_id).height_stack(:,:,1);
mh = nanmax(stack(:));


overlay_fig = figure;%('Visible','off');
hold on
cmap = colormap('pink');
cmap = cmap(75:end,:);
colormap(cmap);

height_slice_1 = imgaussfilt(repo_flux_struct(plot_id).height_stack(:,:,date_ids(1)),1);
surf(height_slice_1)   

% height_slice_2 = imgaussfilt(repo_flux_struct(plot_id).height_stack(:,:,date_ids(2)),1);
% surf(height_slice_2,'FaceAlpha',0.0,'EdgeAlpha',0.5,'EdgeColor',[.75 .75 .75])   
xlabel('x distance (meters)')
ylabel('y distance (meters)')
zlabel('inferred heigh (meters)')
    % h = colorbar;
    % ylabel(h,'plateau height (meters)')

zlim([0 6])
xlim([0 size(height_slice_1,2)])
ylim([0 size(height_slice_1,1)])
view(225,35)
set(gca,'Fontsize',10)
grid on

overlay_fig.Renderer='Painters';
saveas(overlay_fig,[FigDir 'pile_reconstruction_mesh_overlay01.png'])
saveas(overlay_fig,[FigDir 'pile_reconstruction_mesh_overlay01.pdf'])

overlay_fig = figure;%('Visible','off');
hold on
cmap = colormap('pink');
cmap = cmap(75:end,:);
colormap(cmap);

% height_slice_1 = imgaussfilt(repo_flux_struct(plot_id).height_stack(:,:,date_ids(1)),1);
% surf(height_slice_1)   

height_slice_2 = imgaussfilt(repo_flux_struct(plot_id).height_stack(:,:,date_ids(2)),1);
surf(height_slice_2,'FaceAlpha',0.1,'EdgeAlpha',0.5,'EdgeColor',[.75 .75 .75])   
xlabel('x distance (meters)')
ylabel('y distance (meters)')
zlabel('inferred heigh (meters)')
    % h = colorbar;
    % ylabel(h,'plateau height (meters)')

zlim([0 6])
xlim([0 size(height_slice_1,2)])
ylim([0 size(height_slice_1,1)])
view(225,35)
set(gca,'Fontsize',10)
grid on

overlay_fig.Renderer='Painters';
saveas(overlay_fig,[FigDir 'pile_reconstruction_mesh_overlay02.png'])
saveas(overlay_fig,[FigDir 'pile_reconstruction_mesh_overlay02.pdf'])

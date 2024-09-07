clear
close all

DataPath = '../dat/';
FileName = 'all_gangareps_4_23.csv';
RefName = 'Repository_Stats_43examples';
FigDir = '../fig/';
mkdir(FigDir)
% load table
pointTable = readtable([DataPath FileName]);
infoTable = readtable([DataPath RefName]);

%%
rm_filter = infoTable.plateau_area./infoTable.base_area<0.2; % remove cases when very little plateua area was observed
h_v_a = fitlm(infoTable.base_area(~rm_filter), infoTable.max_height(~rm_filter));
coeff_vec = [h_v_a.Coefficients.Estimate' h_v_a.RMSE];

% extract vectors
x_vec = pointTable.POINT_X;
y_vec = pointTable.POINT_Y;
a_vec = pointTable.area;
o_vec = pointTable.ORIG_FID;

% generate unique area vec
region_index = unique(o_vec);
a_vec_u = NaN(size(region_index));
for r = 1:length(region_index)
  a_vec_u(r) = mean(a_vec(o_vec==region_index(r)));
end

% renormalize
norm_factor = 1e4;

min_x = min(x_vec);
min_y = min(y_vec);

max_x = ceil((max(x_vec)-min_x)/norm_factor)*norm_factor;
max_y = ceil((max(y_vec)-min_y)/norm_factor)*norm_factor;

region_index = unique(o_vec);

x_shift = (x_vec - min_x);
y_shift = (y_vec - min_y);

% make blank image
% fov_image = zeros(max_y,max_x);

% add footprints
patch_fig = figure;
hold on
cmap1 = brewermap(length(region_index),'Spectral');
for r = 1:length(region_index)
    region_filter = o_vec==region_index(r);
    patch(x_shift(region_filter),y_shift(region_filter),cmap1(r,:));
end    

area_hist = figure;
histogram(a_vec,'Normalization','probability','FaceColor',cmap1(2,:))

grid on
ylabel('probability')
xlabel('repository surface area (m^2)')
set(gca,'Fontsize',14)

saveas(area_hist,[FigDir 'repo_area_hist.png'])
saveas(area_hist,[FigDir 'repo_area_hist.pdf'])

area_hist = figure;
histogram(a_vec(a_vec<=2e4),'Normalization','probability','FaceColor',cmap1(2,:))

grid on
ylabel('probability')
xlabel('repository surface area (m^2)')
set(gca,'Fontsize',14)
saveas(area_hist,[FigDir 'repo_area_hist_truncated.png'])
saveas(area_hist,[FigDir 'repo_area_hist_truncated.pdf'])
%% Calculate approximate volume for each footprint
close all
% set parameters
angle_of_repose_vec = [20 25 30]; % degrees
m_factor = 1e0;
min_h = min(infoTable.max_height(~rm_filter));
max_h = max(infoTable.max_height(~rm_filter));
n_boots = 50;
region_struct = struct;

% pre-draw heights for bootstraps
for r = 1:length(region_index)    
    region_filter = o_vec==region_index(r);
    for a = 1:length(angle_of_repose_vec)
        % randomly draw height
        max_height_mean = (coeff_vec(1) + coeff_vec(2) * mean(a_vec(region_filter)))*tand(angle_of_repose_vec(a))/tand(30);
        height_dist = makedist('normal', max_height_mean, coeff_vec(3));
        height_dist = truncate(height_dist,min_h,max_h);
        % draw
        region_struct(r).max_height(a,:) = random(height_dist,1,n_boots);
    end
end

% iterate through each region
for r = 1:length(region_index)
  
    % extract points
    region_filter = o_vec==region_index(r);
    x_points = x_shift(region_filter);
    y_points = y_shift(region_filter);            
   
    % readjust scale 
    x_re = m_factor*(x_points-min(x_points));
    y_re = m_factor*(y_points-min(y_points));
    xm = max(x_re);
    ym = max(y_re);

    % generate mask
    cpts = convhull(x_re,y_re);
    BW = poly2mask(1e1 + x_re(cpts), 1e1 + y_re(cpts),2e1+round(ym),2e1+round(xm));

    % calculate distance to edge
    edge_dist_array = bwdist(~BW);
            
    % initialize
    region_struct(r).height_array = NaN(size(BW,1),size(BW,2),length(angle_of_repose_vec));
    region_struct(r).region_vol = NaN(length(angle_of_repose_vec),n_boots);
    
    for a = 1:length(angle_of_repose_vec)        
        for n = 1:n_boots
          
            dist_to_top = region_struct(r).max_height(a,n)/tand(angle_of_repose_vec(a));
          
            height_array = zeros(size(edge_dist_array));
            height_array(edge_dist_array>=dist_to_top) = region_struct(r).max_height(a,n);
            height_array(edge_dist_array<dist_to_top) = edge_dist_array(edge_dist_array<dist_to_top) ...
                          * tand(angle_of_repose_vec(a));

            % turn this into a volume estimate
            region_struct(r).region_vol(a,n) = sum(height_array(:));
            
            if n == 1
                region_struct(r).height_array(:,:,a) = height_array;
            end                        
        end        
    end
    
    region_struct(r).region_area_inf = sum(sum(region_struct(r).height_array(:,:,1)~=0));
    region_struct(r).region_area_actual = mean(a_vec(region_filter));
    region_struct(r).x_points = x_points;
    region_struct(r).y_points = y_points;    
    region_struct(r).angle_of_repose_vec = angle_of_repose_vec;
    region_struct(r).region_id = region_index(r);
end

%%
close all

vol_total_array = NaN(n_boots,length(angle_of_repose_vec));
vol_ind_array = NaN(n_boots,length(region_index),length(angle_of_repose_vec));
for a = 1:length(angle_of_repose_vec)
    total_vol = zeros(1,n_boots);
    for r = 1:length(region_index)
        total_vol = total_vol + region_struct(r).region_vol(a,:);
        vol_ind_array(:,r,a) = region_struct(r).region_vol(a,:);
    end
    vol_total_array(:,a) = total_vol;
end

% total vol
total_vol_fig = figure;

cmap1 = brewermap(length(region_index),'set2');
hold on
for a = 1:length(angle_of_repose_vec)
    scatter(repelem(a,n_boots)+0.25*rand(1,n_boots)-0.125,vol_total_array(:,a),...
      'MarkerfaceColor',cmap1(a+1,:),'MArkerEdgeColor','k','MarkerFaceAlpha',0.05,'MarkerEdgeAlpha',0.05)
    errorbar(a,mean(vol_total_array(:,a)),std(vol_total_array(:,a)),'o','Color','k')
    scatter(a,mean(vol_total_array(:,a)),50,'s',...
      'MarkerfaceColor',cmap1(a+1,:),'MarkerEdgeColor','k','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)
end
ylim([4.5e6 7.75e6])
grid on
ylabel('estimated total storage capacity (m^3)')
xlabel('angle of repose (degrees)')
set(gca,'Fontsize',14)
set(gca,'xtick',1:3,'xticklabel',[20 25 30])

saveas(total_vol_fig,[FigDir 'total_volume_vs_repose_angle.png'])
saveas(total_vol_fig,[FigDir 'total_volume_vs_repose_angle.pdf'])

%% Volume hist
close all

for a = 1:length(angle_of_repose_vec)
  
    total_vol_fig = figure;

    cmap1 = brewermap(length(region_index),'set2');
    hold on
    histogram(vol_ind_array(:,:,a),'Normalization','probability','FaceColor',cmap1(3,:))

    grid on
    xlabel('estimated repository volume (m^3)')
    ylabel('probability')
    set(gca,'Fontsize',14)
%     set(gca,'xtick',1:3,'xticklabel',[20 25 30])
    xlim([0 1e5])
    saveas(total_vol_fig,[FigDir 'volume_hist_rep_angle' num2str(angle_of_repose_vec(a)) '.png'])
    saveas(total_vol_fig,[FigDir 'volume_hist_rep_angle' num2str(angle_of_repose_vec(a))  '.pdf'])
end

%% make sample height array
r_ind = 100;
a_ind = 2;

ha_fig = figure;

cmap2 = flipud(brewermap([],'Blues'));
colormap(cmap2);

imagesc((region_struct(r_ind).height_array(:,:,a_ind)~=0))

xlabel('x distance (meters)')
ylabel('y distance (meters)')

set(gca,'Fontsize',14)

saveas(ha_fig,[FigDir 'footprint_mask.png'])
saveas(ha_fig,[FigDir 'footprint_mask.pdf'])

% example height heatmap
ha_fig = figure;

cmap2 = flipud(brewermap([],'Spectral'));
colormap(cmap2);

imagesc(bwdist(region_struct(r_ind).height_array(:,:,a_ind)==0))

xlabel('x distance (meters)')
ylabel('y distance (meters)')
h = colorbar;
ylabel(h,'distance from pile boundary')
set(gca,'Fontsize',14)

saveas(ha_fig,[FigDir 'edge_distance_calc.png'])
saveas(ha_fig,[FigDir 'edge_distance_calc.pdf'])

% example height heatmap
ha_fig = figure;

cmap2 = flipud(brewermap([],'Spectral'));
colormap(cmap2);

surf(region_struct(r_ind).height_array(:,:,a_ind))

xlabel('x distance (meters)')
ylabel('y distance (meters)')
zlabel('inferred heigh (meters)')
% h = colorbar;
% ylabel(h,'plateau height (meters)')

zlim([0 10])

set(gca,'Fontsize',10)

saveas(ha_fig,[FigDir 'example_pile_reconstruction.png'])
saveas(ha_fig,[FigDir 'example_pile_reconstruction.pdf'])

% scatter(repelem(2,n_boots)+0.5*rand(1,n_boots)-0.25,vol_total_array(:,2))
% scatter(repelem(3,n_boots)+0.5*rand(1,n_boots)-0.25,vol_total_array(:,3))
% close all
% figure;
% imagesc(height_array)
% hold on
% scatter(x_points,y_points,'s')
% % scatter(x_points(BW),y_points(BW),'s')
% scatter(rectx,recty,'d')
% scatter(pgx,pgy,'^')


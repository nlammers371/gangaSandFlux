% Script to make visualization of cluster of repos
clear
close all

DataPath = '../dat/20220214/';
FileName = 'sbapts_02142022putm.csv';
OutDir = '../out/20220214/';
mkdir(OutDir)
FigDir = '../fig/20220214/';
% load table
refTable = readtable([DataPath FileName]);
keyTable = readtable([OutDir 'key_table.csv']);

%%
regionNameCell = {'w_Mirzanagar'};
angle_rep_index = 2;

for i = 1%:length(regionNameCell)
    repo_save_path = [FigDir regionNameCell{i} filesep];
    mkdir(repo_save_path)
    
    % get list of repos falling in this region
    region_repo_filter = contains(keyTable.names_raw,regionNameCell{i});
    region_repo_filter_full = contains(refTable.Name,regionNameCell{i});

    % get list of unique repos and dates
    date_index = unique(keyTable.date_num(region_repo_filter));
    date_repo_array = unique([keyTable.rep_id(region_repo_filter), keyTable.date_num(region_repo_filter)],'rows');
    
    % specify IDs for region to visualize
    regions_to_plot = regionNameCell{i};        

    % load reference set for repo heights
%     load([ProcessedDataPath 'height_fit_struct.mat'],'height_fit_struct')
%     angle_rep_vec = [height_fit_struct.angle_rep];
%     angle_of_repose = angle_rep_vec(angle_rep_index);

    % extract vectors
    x_vec_ft = refTable.POINT_X(region_repo_filter_full);
    y_vec_ft = refTable.POINT_Y(region_repo_filter_full);   
    names_ft = refTable.Name(region_repo_filter_full);   
    
    % renormalize position data
    norm_factor = 1e1;

    min_x = floor(min(x_vec_ft)/norm_factor)*norm_factor;
    min_y = floor(min(y_vec_ft)/norm_factor)*norm_factor;

    max_x = ceil((max(x_vec_ft)-min_x)/norm_factor)*norm_factor;
    max_y = ceil((max(y_vec_ft)-min_y)/norm_factor)*norm_factor;

    x_shifted = (x_vec_ft - min_x);
    y_shifted = (y_vec_ft - min_y);

    % initialize empty frame to populate
    region_fp_array = zeros(max_y,max_x,length(date_index));

    rng(102);
    for d = 1:length(date_index)
      
        fp_fig = figure('Visible','off');
        cmap = colormap('pink');
        
        hold on
        
        region_fp_frame = zeros(max_y,max_x);
        regions_to_plot = date_repo_array(date_repo_array(:,2)==date_index(d),1);
        
        for r = 1:length(regions_to_plot)
                   
            repo_name = keyTable.names_raw(keyTable.rep_id==regions_to_plot(r));
            name_filter = ismember(names_ft,repo_name(1));
            
            % extract footprint and calculate total volume    
            xfp = x_shifted(name_filter);
            yfp = y_shifted(name_filter);
            
            patch(xfp, yfp, cmap(150,:),'EdgeColor','k')
            
            % generate mast
%             repo_mask = poly2mask(xfp, yfp, max_y, max_x);    
%             region_fp_frame = region_fp_frame | repo_mask;           
            
            
%                                  
        end
        region_fp_array(:,:,d) = region_fp_frame;   
        region_fp_frame3 = double(1*region_fp_frame);
        dist_frame = bwdist(~region_fp_frame3);
        region_fp_frame3(region_fp_frame==1) = 2;
        region_fp_frame3(region_fp_frame==0) = 1;
        region_fp_frame3(round(dist_frame)==1) = 0;
                        
%         cmap = colormap('pink');
%         cmap = cmap([1 50 150],:);
%         colormap(cmap);
%         
%         p = imagesc(flipud(region_fp_frame3));
%         set(p, 'EdgeAlpha', 0.1);

        xlabel('x position (meters)')
        ylabel('y position (meters)')
        box on
        set(gca,'FontSize',14)
        fp_fig.Renderer='Painters';
        
        saveas(fp_fig,[repo_save_path 'ROI_reconstruction_' datestr(datetime(date_index(d),'ConvertFrom','datenum')) '.png'])
        saveas(fp_fig,[repo_save_path 'ROI_reconstruction_' datestr(datetime(date_index(d),'ConvertFrom','datenum')) '.pdf'])
        
    end
    
end



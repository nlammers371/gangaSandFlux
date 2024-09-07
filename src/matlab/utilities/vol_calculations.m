function [height_array, dist_array, repo_mask, x_bounds, y_bounds, max_height_vals] = ...
        vol_calculations(x_vec,y_vec,fit_struct,x_bounds,y_bounds,SA,n_boots, max_height_vals)    
                  
    % extract values
    angle_of_repose = fit_struct.angle_rep;    
      
    % check for pre assigned image bounds
    if isempty(x_bounds)
        x_bounds(1) = min(x_vec);
        x_bounds(2) = max(x_vec);
        y_bounds(1) = min(y_vec);
        y_bounds(2) = max(y_vec);        
    end
    
    % readjust scale 
    xre = (x_vec-x_bounds(1));
    yre = (y_vec-y_bounds(1));
    
    xm = x_bounds(2)-x_bounds(1);
    ym = y_bounds(2)-y_bounds(1);
    
    % generate region mask     
    repo_mask = poly2mask(1e1 + xre, 1e1 + yre,2e1+round(ym),2e1+round(xm));
    if isempty(SA)
        SA = sum(repo_mask(:));
    end
    % generate distance map
    dist_array = bwdist(~repo_mask);
    
    % calculate max height (assume deterministic for now) 
    if isempty(max_height_vals)
        mh_mean = fit_struct.x_fit(1) * SA ./ (fit_struct.x_fit(2) + SA);
        mh_sigma = fit_struct.mdl_se;
        min_height = min(fit_struct.max_height_vec);
        height_dist = makedist('normal', mh_mean, mh_sigma);
        height_dist = truncate(height_dist,min_height,Inf);    
        max_height_vals = random(height_dist,1,n_boots);
    end
    
    % calculate repo capacity   
    height_array = zeros(size(dist_array,1),size(dist_array,2),n_boots);
    for n = 1:n_boots
        dist_to_top = max_height_vals(n)/tand(angle_of_repose);
        slice = height_array(:,:,n);
        slice(dist_array>=dist_to_top) = max_height_vals(n);
        slice(dist_array<dist_to_top) = dist_array(dist_array<dist_to_top) ...
                                                    * tand(angle_of_repose);
        height_array(:,:,n) = slice;
    end
    
% Preprocess - coordinating

% coordinates for design variable nodes
% rectangle
tp=0:par.hx:par.hx*par.nelx; tp = repmat(tp, par.nely+1, 1); dv_coords=tp(:);
tp=par.hy*par.nely:-par.hy:0; tp = repmat(tp, 1, par.nelx+1); dv_coords=[dv_coords, tp(:)];

% outer fixed points
mins = min(dv_coords, [], 1); maxs = max(dv_coords, [], 1);
min_x = mins(1); min_y = mins(2); max_x = maxs(1); max_y = maxs(2);
nx_cover = par.nelx+1; ny_cover = par.nely+1;
cover_coords = [];
for i=1:par.coverthickness    
    cover_one = [ (min_x - par.hx)*ones(ny_cover, 1), (max_y:-par.hy:min_y)';
        min_x-par.hx, min_y-par.hy;
        (min_x:par.hx:max_x)', (min_y-par.hy)*ones(nx_cover, 1);
        max_x+par.hx, min_y-par.hy;
        max_x+par.hx*ones(ny_cover, 1), (min_y:par.hy:max_y)';
        max_x+par.hx, max_y+par.hy;
        (max_x:-par.hx:min_x)', (max_y+par.hy)*ones(nx_cover, 1);
        min_x-par.hx, max_y+par.hy];
    cover_coords = [cover_coords; cover_one];
    mins = mins-[par.hx, par.hy]; maxs = maxs+[par.hx, par.hy];
    min_x = mins(1); min_y = mins(2); max_x = maxs(1); max_y = maxs(2);
    nx_cover = nx_cover+2; ny_cover = ny_cover+2;
end
len_cover = size(cover_coords,1);

disp('Coordinates computed');
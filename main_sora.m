
files = dbstack;
curpath = fullfile(pwd,'results',[datestr(now, 'yymmdd_HH;MM;SS'), files(end).name]);
mkdir(fullfile(curpath, 'figures'))
prob = 1;

%% Problem definition
if prob == 1    
    % Domain size
    nelx = 160; nely = 40;
    hx = 1; hy = 1;
    % Load definition
    load_coords = [nelx, nely/2];
    load_vector = [0, -1];
    % Boundary condition
    fixed_coords = [0, nan, 1; 0, nan, 2]; % nan: all
    % Constraint
    c_target = 500;
    beta_target = 1;
    dtostep = 1;        
elseif prob == 2    
    % % Domain size
    nelx = 100; nely = 50;
    hx = 1; hy = 1;
    % Load definition
    load_coords = [0, 0];
    load_vector = [1, 0];
    % Boundary condition
    fxl = [0:0.5:3]'; mfxl = max(fxl);
    ll=length(fxl);
    fixed_coords = [zeros(ll, 1), nely-fxl, ones(ll, 1); zeros(ll, 1), nely-fxl, 2*ones(ll, 1); ...
        nan, 0, 2]; % nan: all
    % For compliant mechanism
    output_height = 15;
    compliant_output_coords = [nelx, output_height, 2]; % x, y, xy-dir, direction (+, -)
    spring_coords = [0, 0, 1, 1; nelx, output_height, 2, 0.01]; % x, y, xy-dir, spring coefficient
    % Constraint
    c_target = -0.4;
    beta_target = 1;
    dtostep = 1;    
end

% Uncertainty
sigma_corr = ones((nelx+1)*(nely+1), 1);
error_target = 0.01;
corr_length = 80;
% etc.
dtomode = 0;
candi_additional_dist = 1;
beta_end = 32;

%% Preprocesses
preprocess_parameter
preprocess_coordinating
preprocess_fem
preprocess_candi
preprocess_kl

if prob == 1
    if ~isempty(cover_coords)
        cover_coords = cover_coords(find(cover_coords(:, 1)>=0), :);
    end
    len_cover = length(cover_coords);
    par.len_cover1 = len_cover; len_cover2 = 0;
    density_coords = [dv_coords; cover_coords];
    rho_cover = zeros(len_cover, 1);
    rho_dv = ones(par.dv_num, 1);
    dv_ind = [1:par.dv_num]';
    
elseif prob == 2    
    cover_ind =((cover_coords(:, 1)<0)+(cover_coords(:, 1)>nelx)+(cover_coords(:, 2)>0));
    cover_coords1 = cover_coords(find(cover_ind), :);
    cover_coords2 = cover_coords(find(~cover_ind), :);    
    cover_coords = [cover_coords1; cover_coords2];
    len_cover1 = size(cover_coords1, 1);
    len_cover2 = size(cover_coords2, 1);    
    len_cover = len_cover1 + len_cover2;
    par.len_cover1 = len_cover1;
    density_coords = [dv_coords; cover_coords];
    
    rho_cover = zeros(len_cover, 1);     
    rho_cover(find( (cover_coords(:, 1)<0) .* ...
        ( (cover_coords(:, 2)>=nely-mfxl) .* (cover_coords(:, 2)<=nely)))) = 1;
    
    central_fix = 2;
    rho_cover(find((cover_coords(:, 2)<=central_fix).*(cover_coords(:, 1)<0))) = 1;    
    dv_ind = find((dv_coords(:, 1)<nelx-output_height*2)+(dv_coords(:, 2)>output_height));
    rho_dv = ones(par.dv_num, 1);
    rho_dv(find((dv_coords(:, 2)==output_height).*(dv_coords(:, 1)>=nelx-output_height*2)))=1;
    rho_dv(find((dv_coords(:, 2)<output_height).*(dv_coords(:, 1)>=nelx-output_height*2))) =par.rhomin;
end

preprocess_findsx
preprocess_filter
move=0.2;

%% Optimize

[rho_final, time_opt, total_iter] = ...
    optimize_sora(rho_dv, rho_cover, beta_target, c_target, par, dv_coords, cover_coords, len_cover, density_coords, ...
    fem, xs_quads, sx_cover, sx_default, sx_pseudograd, kl, candi_dv, candi_num, dtomode, dv_ind, H, Hs, curpath, move);




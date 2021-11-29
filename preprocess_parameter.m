% Preprocess - parameter


% Dimensional parameters
par.dim = 2;
par.hx = hx; par.hy = hy;
par.nelx = nelx; par.nely = nely;
par.lx = par.nelx*par.hx; par.ly = par.nely*par.hy;
par.coverthickness = 2;

% Uncertainty parameters
par.corr_length = 80;
par.error_target_KL = 0.1;
par.eignum = 20; % initial guess for the dimension of discretized random field
par.dtostep = dtostep-1;

% Material parameters
par.E = 1; par.nu = 0.3; 

D0 = par.E/(1-par.nu^2) * [1 par.nu 0; par.nu 1 0; 0 0 0.5*(1-par.nu)]; % plane stress
% D0 = E/((1+par.nu)*(1-2*par.nu)) * [1-par.nu par.nu 0; par.nu 1-par.nu 0; 0 0 (1-2*par.nu)/2]; % plane strain


% Optimization parameters
par.dv_num = (nelx+1)*(nely+1);
par.rhomin = 1e-3;
par.penal = 3;

par.tol_do = 2e-2;
par.tol_ra = 2e-2;
par.Rc = 2;
par.Rc_pseudograd = 2;
par.l_candi = par.Rc+candi_additional_dist;
par.dvmax_candi =candi_additional_dist;




switch lower(interpolation)
    case 'shepard'
        par.intp_func = 1;
    case 'exp_1'
        par.intp_func = 2;
    case 'exp_2'
        par.intp_func = 3;
    case 'shepard_1'
        par.intp_func = 4;
    case 'shepard_3'
        par.intp_func = 5;
    case 'shepard_5'
        par.intp_func = 6;
    case 'exp_3'
        par.intp_func = 7;
    case 'exp_inv_1'
        par.intp_func = 8;
    case 'exp_inv_2'
        par.intp_func = 9;
    case 'exp_inv_3'
        par.intp_func = 10;
    case 'log_inv_1'
        par.intp_func = 11;
    case 'log_inv_2'
        par.intp_func = 12;
    case 'log_inv_3'
        par.intp_func = 13;
    case 'exp_inv_0_2'
        par.intp_func = 14;
    case 'exp_inv_0_1'
        par.intp_func = 15;
    case 'csch'
        par.intp_func = 16;
    case 'coth'
        par.intp_func = 17;
        
    otherwise
        error('Check interpolation')
end


[grid_x, grid_y] = meshgrid(-par.Rc:par.hx:par.Rc, -par.Rc:par.hy:par.Rc);
grid_in_rc = [grid_x(:), grid_y(:)];
dist_in_rc = vecnorm(grid_in_rc, 2, 2);
inrc = intersect(find(dist_in_rc>0), find(dist_in_rc<par.Rc));
grid_in_rc = grid_in_rc(inrc, :);
dist_in_rc = dist_in_rc(inrc, :);

par.Mp = sum(abs(grid_in_rc(:,1)) .* abs(interp_function(dist_in_rc, par.intp_func, 1)) ./dist_in_rc) ...
    / 2 / sum(interp_function(dist_in_rc, par.intp_func, 0));




disp('Parameter defined');
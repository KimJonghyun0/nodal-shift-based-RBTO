% Preprocess - parameter
par.prob = prob;
% Dimensional parameters
par.dim = 2;
par.hx = hx; par.hy = hy;
par.nelx = nelx; par.nely = nely;
par.lx = par.nelx*par.hx; par.ly = par.nely*par.hy;
par.coverthickness = 2;
% Uncertainty parameters
par.corr_length = corr_length;
par.eignum = 20;
par.dtostep = dtostep-1;
% Optimization parameters
par.dv_num = (nelx+1)*(nely+1);
par.rhomin = 0;
par.penal = 3;
par.beta_end = beta_end;
par.tol_do = 1e-2;
par.tol_ra = 1e-2;
par.Rc = 2;
par.Rc_pseudograd = 2;
par.l_candi = par.Rc+candi_additional_dist;
par.dvmax_candi =candi_additional_dist;

% Material parameters
par.E = 1; par.nu = 0.3;
par.Emin = par.E*1e-3^par.penal;
D0 = par.E/(1-par.nu^2) * [1 par.nu 0; par.nu 1 0; 0 0 0.5*(1-par.nu)]; % plane stress
[grid_x, grid_y] = meshgrid(-par.Rc:par.hx:par.Rc, -par.Rc:par.hy:par.Rc);
grid_in_rc = [grid_x(:), grid_y(:)];
distsq_in_rc = sum(grid_in_rc.^2, 2);
inrc = intersect(find(distsq_in_rc>0), find(distsq_in_rc<par.Rc.^2));
grid_in_rc = grid_in_rc(inrc, :);
distsq_in_rc = distsq_in_rc(inrc, :);
par.Mp = sum(abs(grid_in_rc(:,1)) .* distsq_in_rc.^(-2)) ...
    / sum(1./distsq_in_rc);
disp('Parameter defined');
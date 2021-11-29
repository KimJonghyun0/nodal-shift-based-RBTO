function [sx_mpp, xj_moved_mpp, xi_mpp]=reliability_analysis(rho, beta_target, par, fem, kl, pseudograd_rho, dv_coords, sx_cover, xs_quads, candi_dv, candi_num, search_method)
% Reliability analysis

xi = zeros(kl.m, 1); % initial xi

ndcdxii_old = zeros(kl.m, 1);
ndcdxii_cur = zeros(kl.m, 1); % For HMV

crit = 1; % criterion

mpp_iter = 1;

try kl.mean; catch kl.mean = zeros(par.dv_num, 1); end
% Find MPP
while crit>par.tol_ra/10
    yhat = kl.phi * xi + kl.mean;
    dxj = pseudograd_rho .* yhat / par.Mp;
    xj_moved = dv_coords + dxj;
    
    % Tree search or Direct search
    if strcmp(search_method, 'tree')
        sx_xi = findsx_tree(xj_moved, tree, par.Rc, fem.numelem, xs_quads_matrix, sx_cover);
    elseif strcmp(search_method, 'direct')
        sx_xi = findsx_direct(xs_quads, xj_moved, par.Rc, sx_cover);
    elseif strcmp(search_method, 'candi')
        sx_xi = findsx_candi(xs_quads, xj_moved, par.Rc, par.dvmax_candi, dxj, candi_dv, candi_num, par.intp_func, sx_cover);
    end
    % FEM
    U_xi = structural_analysis(fem, par, sx_xi, rho);
    % Senstivity
    dcdxii = zeros(1, kl.m);
    for elem = 1:fem.numelem
        ue = U_xi(fem.edof(elem, :), 1); % KU=F
        ve = U_xi(fem.edof(elem, :), end); % KV=L, for compliant mechanism
        penal_bdbwa = par.penal*fem.bdbwa(:, :, elem);
        for quad_pt = 1:3
            rho_qd = sx_xi.interp(1:sx_xi.num(elem, quad_pt), elem, quad_pt)'...
                *rho(sx_xi.nodes_in_sx(1:sx_xi.num(elem, quad_pt), elem, quad_pt));
            under_dvnum = (sx_xi.nodes_in_sx(1:sx_xi.num(elem, quad_pt), elem, quad_pt)<=par.dv_num);
            nodes_in_sx_domain = sx_xi.nodes_in_sx(under_dvnum, elem, quad_pt);
            xj_subs_xg_domain = -sx_xi.del_xs(under_dvnum, :, elem, quad_pt);          
            dist_domain = sx_xi.dist(under_dvnum, elem, quad_pt);
            df = interp_function(dist_domain, par.intp_func, 1);
            drho_dxii_xj = sum(pseudograd_rho(nodes_in_sx_domain, :).*xj_subs_xg_domain,2).*df./dist_domain ...
                .*(rho(nodes_in_sx_domain)-rho_qd)./sx_xi.denominator(elem, quad_pt) .* kl.phi(nodes_in_sx_domain, :);
            drho_dxii = sum(drho_dxii_xj, 1);
            dcdxii = dcdxii - drho_dxii*(ve'*rho_qd^(par.penal-1) * penal_bdbwa*ue);
        end
        if isnan(dcdxii)
            error('NaN dcdxii')
        end
    end
    dcdxii = dcdxii / par.Mp; % normalization
    
    % Update
    ndcdxii_old2 = ndcdxii_old; ndcdxii_old = ndcdxii_cur ;
    dcdxii_cur = dcdxii(:);
    ndcdxii_cur  = dcdxii_cur(:)./norm(dcdxii_cur);
    sigma_hmv = (ndcdxii_cur  - ndcdxii_old)' * (ndcdxii_old - ndcdxii_old2);
    if sigma_hmv > 0 || mpp_iter<3
        alpha = ndcdxii_cur ;
    else
        alpha = (ndcdxii_cur +ndcdxii_old+ndcdxii_old2)/norm(ndcdxii_cur +ndcdxii_old+ndcdxii_old2);
    end
    xi_new = beta_target * alpha;
    xi_change = max(abs(xi-xi_new));
    crit = min(abs(alpha' * xi/norm(xi) - 1), 1);
    xi=xi_new;
    mpp_iter = mpp_iter + 1;
    
    if isnan(xi)
        error('NaN xi')
    end
    if mpp_iter>100
        warning('Too many mpp searching iterations')
        break;
    end
end


xi_mpp = xi;
yhat_mpp = kl.phi * xi_mpp;
dxj_mpp = pseudograd_rho .* yhat_mpp / par.Mp;
xj_moved_mpp =  dv_coords + dxj_mpp;


if strcmp(search_method, 'tree')
    sx_mpp = findsx_tree(xj_moved_mpp, tree, par.Rc, fem.numelem, xs_quads_matrix, sx_cover);
elseif strcmp(search_method, 'direct')
    sx_mpp = findsx_direct(xs_quads, xj_moved_mpp, par.Rc, sx_cover);
elseif strcmp(search_method, 'candi')
    sx_mpp = findsx_candi(xs_quads, xj_moved, par.Rc, par.dvmax_candi, dxj, candi_dv, candi_num, par.intp_func, sx_cover);
end

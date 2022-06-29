function [sx_xi, xj_moved, xi, U_xi, mpp_iter]=reliability_analysis(rho, beta_target, par, fem, kl, pseudograd_rho, dv_coords, sx_cover, xs_quads, candi_dv, candi_num)
xi_new = zeros(kl.cumsum_m(end), 1);
ndcdxii_old = zeros(kl.cumsum_m(end), 1);
ndcdxii_cur = zeros(kl.cumsum_m(end), 1);
crit = 1;
mpp_iter = 0;
yhat = zeros(par.dv_num, kl.dim);
try kl.mean; catch kl.mean = zeros(par.dv_num, 1); end
% Find MPP
while crit>par.tol_ra/10
    mpp_iter = mpp_iter + 1;
    xi=xi_new;    
    for i_dim = 1:kl.dim
        yhat(:, i_dim) = kl.phi(:, kl.cumsum_m(i_dim)+1:kl.cumsum_m(i_dim+1)) * xi(kl.cumsum_m(i_dim)+1:kl.cumsum_m(i_dim+1));
    end
    yhat = yhat + kl.mean;    
    dxj = pseudograd_rho .* yhat / par.Mp;
    xj_moved = dv_coords + dxj;
    sx_xi = findsx_candi(xs_quads, xj_moved, par.Rc, par.dvmax_candi, dxj, candi_dv, candi_num, sx_cover);
    % FEM
    U_xi = structural_analysis(fem, par, sx_xi, rho);
    % Senstivity
    dcdxii = zeros(1, kl.cumsum_m(end));
    for elem = 1:fem.numelem
        ue = U_xi(fem.edof(elem, :), 1); % KU=F
        ve = U_xi(fem.edof(elem, :), end); % KV=L, for compliant mechanism
        dE_penal_bdbwa = (par.E-par.Emin)*par.penal*fem.bdbwa(:, :, elem);
        for quad_pt = 1:3
            rho_qd = sx_xi.interp(1:sx_xi.num(elem, quad_pt), elem, quad_pt)'...
                *rho(sx_xi.nodes_in_sx(1:sx_xi.num(elem, quad_pt), elem, quad_pt));
            under_dvnum = (sx_xi.nodes_in_sx(1:sx_xi.num(elem, quad_pt), elem, quad_pt)<=par.dv_num);
            nodes_in_sx_domain = sx_xi.nodes_in_sx(under_dvnum, elem, quad_pt);
            zg_subs_xj_domain = sx_xi.del_xs(under_dvnum, :, elem, quad_pt);
            dist_sq_domain = sx_xi.dist_sq(under_dvnum, elem, quad_pt);
            drho_dxii_xj = 2*sum(pseudograd_rho(nodes_in_sx_domain, :).*zg_subs_xj_domain,2)./dist_sq_domain.^2 ...
                .*(rho(nodes_in_sx_domain)-rho_qd)./sx_xi.denominator(elem, quad_pt) .* kl.phi(nodes_in_sx_domain, :);
            drho_dxii = sum(drho_dxii_xj, 1);
            dcdxii = dcdxii - drho_dxii*(ve'*rho_qd^(par.penal-1) * dE_penal_bdbwa*ue);
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
    crit = min(abs(alpha' * xi/norm(xi) - 1), 1);    
end
end

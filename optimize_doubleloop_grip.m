function [rho_dv]=optimize_doubleloop_grip(rho_dv, rho_cover, beta_target, c_target, par, dv_coords, cover_coords, len_cover, density_coords, ...
    fem, xs_quads, sx_cover, sx_default, sx_pseudograd, kl, candi_dv, candi_num, search_method, dtomode, dv_ind)


[m_MMA, n_MMA, iter, xval, xold1, xold2, low, upp, xmin, xmax, a0_MMA, a_MMA, c_MMA, d_MMA, MMA_movelimit] =...
    initialize_MMA(length(dv_ind), par.rhomin);
rho = [rho_dv; rho_cover];

rho(par.dv_num+par.len_cover1+1:end) = rho_dv(sx_pseudograd.ind_sym);  

x_change = 1;
time_start_iter = tic;
c_mpp = 1;
while x_change > par.tol_do && iter < 500 || iter < 10
    pseudograd_rho = pseudograd(par, rho, sx_pseudograd);
    
    % Reliability analysis: Search MPP and perform structural analysis
    if iter > par.dtostep && dtomode == 0 && c_mpp < -0.01
        [sx_mpp, dv_coords_mpp, xi_mpp]=reliability_analysis(rho, beta_target, par, fem, kl, pseudograd_rho, dv_coords, sx_cover, xs_quads, candi_dv, candi_num, search_method);
    else
        sx_mpp = sx_default;
        dv_coords_mpp = dv_coords;
    end
    U_mpp = structural_analysis(fem, par, sx_mpp, rho);
    c_mpp = fem.F(:, end)'*U_mpp(:, 1);


    % Sensitivity analysis
    dcdrho = zeros(1,par.dv_num + len_cover);
    vol_det = 0;
    dvdrho = zeros(1,par.dv_num + len_cover);
    
    for elem = 1:fem.numelem
        ue = U_mpp(fem.edof(elem, :), 1);   
        ve = U_mpp(fem.edof(elem, :), end);   
        for quad_pt = 1:3
            num_mpp_i = sx_mpp.num(elem, quad_pt);
            node_in_elem_qd_mpp = sx_mpp.nodes_in_sx(1:num_mpp_i, elem, quad_pt);
            dcdrho(node_in_elem_qd_mpp) = dcdrho(node_in_elem_qd_mpp) ...
                - par.penal*(sx_mpp.interp(1:num_mpp_i, elem, quad_pt)'*rho(node_in_elem_qd_mpp))^(par.penal-1) ...
                * sx_mpp.interp(1:num_mpp_i, elem, quad_pt)'*(ve'*fem.bdbwa(:, :, elem)*ue);
            
            num_default_i = sx_default.num(elem, quad_pt);
            node_in_elem_qd_default = sx_default.nodes_in_sx(1:num_default_i, elem, quad_pt);
            dvdrho(node_in_elem_qd_default) = dvdrho(node_in_elem_qd_default) ...
                + fem.area(elem)*fem.weights_qd*sx_default.interp(1:num_default_i, elem, quad_pt)';
            vol_det = vol_det + fem.area(elem)*fem.weights_qd ...
                * sx_default.interp(1:num_default_i, elem, quad_pt)'*rho(node_in_elem_qd_default);
        end
    end
    dcdrho = dcdrho(1:par.dv_num);
    dvdrho = dvdrho(1:par.dv_num);
    
    
    % Design update: perform MMA
    
    xval = rho_dv(:); f0val = vol_det; df0dx = dvdrho; fval = c_mpp - c_target; dfdx = dcdrho;
     
    xval = xval(dv_ind); df0dx = df0dx(dv_ind); dfdx = dfdx(dv_ind); 
    iter = iter + 1;
    [xmma,ymma,zmma,lam,xsi,eta_mma,mu,zet,s,low,upp] = ...
        mmasub(m_MMA, n_MMA, iter, xval, xmin, xmax, xold1, xold2, f0val,df0dx(:),fval,dfdx(:)',low,upp,a0_MMA,a_MMA,c_MMA,d_MMA);
    rho_dv(dv_ind) = xmma; xold2=xold1; xold1=xval; rho = [rho_dv; rho_cover];
    rho(par.dv_num+par.len_cover1+1:end) = rho_dv(sx_pseudograd.ind_sym);  
    x_change = max(abs(xmma-xval));
    
    fprintf('Iter.: %d,  Comp.: %7.3f, Vol.: %7.3f, Ch.: %6.3f\n', iter, c_mpp, vol_det, x_change)

    density_coords_mpp = [dv_coords_mpp; cover_coords];

    subplot(2,1,1)
    scatter(density_coords(:, 1), density_coords(:, 2), 5,-rho, 'filled');
    axis equal; axis off; axis tight;
	colormap(winter); grid on;
	caxis([-1,0]); drawnow;
    subplot(2,1,2)
    scatter(density_coords_mpp(:, 1), density_coords_mpp(:, 2),5, -rho, 'filled');
    axis equal; axis off; axis tight;
	caxis([-1,0]); drawnow;

end
    
    
    
    
    
    
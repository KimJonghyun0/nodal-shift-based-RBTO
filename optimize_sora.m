function [rho_dv, time_opt, iter]=optimize_sora(rho_dv, rho_cover, beta_target, c_target, par, dv_coords, cover_coords, len_cover, density_coords, ...
    fem, xs_quads, sx_cover, sx_default, sx_pseudograd, kl, candi_dv, candi_num, dtomode, dv_ind, H, Hs, curpath, move)


[m_MMA, n_MMA, iter, xval, xold1, xold2, low, upp, xmin, xmax, a0_MMA, a_MMA, c_MMA, d_MMA, MMA_movelimit] =...
    initialize_MMA(length(dv_ind), par.rhomin);
rho = [rho_dv; rho_cover];
rho(par.dv_num+par.len_cover1+1:end) = rho_dv(sx_pseudograd.ind_sym);


time_start_iter = tic;
c_mpp = 1;
beta = 1; eta = 0.5;
loopbeta = 0;
tol_in1 = 0.1;
inner_max = 20;
iter_out = 0; iter_total = 0;
yhat = zeros(par.dv_num, kl.dim);
vol_det = 0;
sxdef = 1;
c_change = m_MMA; c_changed = -1;
x_change_out = 1; f_change_out = 1;
beta_changed_or_first_iter = 0;

rho_filtered = (H*rho(:))./Hs;
rho_Phys = ( tanh(beta*eta)+tanh(beta*(rho_filtered(:)-eta)) ) / ( tanh(beta*eta)+tanh(beta*(1-eta)) );
dx_Phys = beta*sech(beta*(rho_filtered(:)-eta)).^2 / ( tanh(beta*eta)+tanh(beta*(1-eta)) );
pseudograd_rho = pseudograd(par, rho_Phys, sx_pseudograd);

while (loopbeta<5 || x_change_out>1e-2 && f_change_out>1e-3) && iter_total < 500 || iter_total < 10 || beta < par.beta_end
    
    %% DO
    [m_MMA, n_MMA, iter, ~, ~, ~, ~, ~, xmin, xmax, a0_MMA, a_MMA, c_MMA, d_MMA, MMA_movelimit] =...
        initialize_MMA(length(dv_ind), par.rhomin);
    
    x_old_out = rho(dv_ind);
    f_old_out = vol_det;
    x_change_in = 1;
    
    while ( loopbeta<5 || x_change_in > tol_in1 || iter < 4 ) && iter < inner_max && c_change > 0 ...
            && iter_total < 500 &&  ~beta_changed_or_first_iter
        c_mpp_old = c_mpp;
        
        dxj = pseudograd_rho .* yhat / par.Mp;
        xj_moved = dv_coords + dxj;
        if sxdef == 0
            sx_mpp = findsx_candi(xs_quads, xj_moved, par.Rc, par.dvmax_candi, dxj, candi_dv, candi_num, sx_cover);
        elseif sxdef == 1
            sx_mpp = sx_default;
        end
        density_coords_mpp = [xj_moved; cover_coords];
        
        U_mpp = structural_analysis(fem, par, sx_mpp, rho_Phys);
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
                % rho_qd = sx_mpp.interp{elem, quad_pt}'*rho(node_in_elem_qd);
                dcdrho(node_in_elem_qd_mpp) = dcdrho(node_in_elem_qd_mpp) ...
                    - par.penal*(par.E-par.Emin)*(sx_mpp.interp(1:num_mpp_i, elem, quad_pt)'*rho_Phys(node_in_elem_qd_mpp))^(par.penal-1) ...
                    * sx_mpp.interp(1:num_mpp_i, elem, quad_pt)'*(ve'*fem.bdbwa(:, :, elem)*ue);
                
                num_default_i = sx_default.num(elem, quad_pt);
                node_in_elem_qd_default = sx_default.nodes_in_sx(1:num_default_i, elem, quad_pt);
                dvdrho(node_in_elem_qd_default) = dvdrho(node_in_elem_qd_default) ...
                    + fem.area(elem)*fem.weights_qd*sx_default.interp(1:num_default_i, elem, quad_pt)';
                vol_det = vol_det + fem.area(elem)*fem.weights_qd ...
                    * sx_default.interp(1:num_default_i, elem, quad_pt)'*rho_Phys(node_in_elem_qd_default);
            end
        end
        
        dcdrho(:) = H*(dcdrho(:).*dx_Phys(:)./Hs);
        dvdrho(:) = H*(dvdrho(:).*dx_Phys(:)./Hs);
        
        dcdrho = dcdrho(1:par.dv_num);
        dvdrho = dvdrho(1:par.dv_num);
        
        
        % Design update: perform MMA
        
        xval = rho_dv(dv_ind); f0val = vol_det; df0dx = dvdrho(dv_ind); fval = c_mpp - c_target; dfdx = dcdrho(dv_ind);
        xmin = max(par.rhomin, xval - move); xmax = min(1, xval + move);
        iter = iter + 1;
        loopbeta = loopbeta + 1;
        iter_total = iter_total+1;
        [xmma,ymma,zmma,lam,xsi,eta_mma,mu,zet,s,low,upp] = ...
            mmasub(m_MMA, n_MMA, iter, xval, xmin, xmax, xold1, xold2, f0val,df0dx(:),fval,dfdx(:)',low,upp,a0_MMA,a_MMA,c_MMA,d_MMA);
        rho_dv(dv_ind) = xmma; xold2=xold1; xold1=xval; rho = [rho_dv; rho_cover];
        rho(par.dv_num+par.len_cover1+1:end) = rho_dv(sx_pseudograd.ind_sym);
        x_change_in = max(abs(xmma-xval));
        
        fprintf('Total iter: %3d, Iter.: %3d,  Comp.: %7.3f, Vol.: %7.3f, Ch.: %6.3f\n', ...
            iter_total, iter, c_mpp, vol_det, x_change_in)
        subplot(2,1,1)
        scatter(density_coords(:, 1), density_coords(:, 2), 5,-rho_Phys, 'filled');
        axis equal; axis off; axis tight;
        colormap(winter); grid on;
        caxis([-1,0]); drawnow;
        subplot(2,1,2)
        scatter(density_coords_mpp(:, 1), density_coords_mpp(:, 2),5, -rho_Phys, 'filled');
        axis equal; axis off; axis tight;
        saveas(gcf, fullfile(curpath, 'figures',['iter',num2str(iter_total),'.png']))
        
        if loopbeta >= 5 && beta < par.beta_end && (loopbeta >= 50 || x_change_in <= par.tol_do) && iter > 3
            beta = min(2*beta, par.beta_end);
            loopbeta = 0;
            x_change_in = 1;
            beta_changed_or_first_iter = 1;
            fprintf('Parameter beta increased to %g.\n',beta);
        end
        
        rho_filtered = (H*rho(:))./Hs;
        rho_Phys = ( tanh(beta*eta)+tanh(beta*(rho_filtered(:)-eta)) ) / ( tanh(beta*eta)+tanh(beta*(1-eta)) );
        dx_Phys = beta*sech(beta*(rho_filtered(:)-eta)).^2 / ( tanh(beta*eta)+tanh(beta*(1-eta)) );
        pseudograd_rho = pseudograd(par, rho_Phys, sx_pseudograd);
        
        if c_changed < 0
            c_change = c_change - sum((c_mpp_old > -0.01) .* (c_mpp < -0.01));
        end
        if iter_total < 2 && par.prob == 1
            beta_changed_or_first_iter = 1;
        end
        
    end
    beta_changed_or_first_iter = 0;
    if c_change == 0, c_changed = 1; c_change = 100; end
    
    
    % Reliability analysis: Search MPP and perform structural analysis
    if iter_out > par.dtostep && dtomode == 0 && (size(fem.F, 2) == 1 || c_mpp < -0.01)
        [~, ~, xi_mpp, ~, mpp_iter]=reliability_analysis(rho_Phys, beta_target, par, fem, kl, pseudograd_rho, dv_coords, sx_cover, xs_quads, candi_dv, candi_num);
        for i_dim = 1:kl.dim
            yhat(:, i_dim) = kl.phi(:, kl.cumsum_m(i_dim)+1:kl.cumsum_m(i_dim+1)) * xi_mpp(kl.cumsum_m(i_dim)+1:kl.cumsum_m(i_dim+1));
        end
        yhat = yhat + kl.mean;
        sxdef = 0;
    else
        yhat = zeros(par.dv_num, kl.dim, m_MMA);
        yhat = yhat + kl.mean;
        sxdef = 1;
        mpp_iter = 0;
    end
    
    iter_out = iter_out + 1;
    
    if tol_in1>par.tol_do
        tol_in1 = max(par.tol_do, tol_in1-0.2);
    end
    
    x_new_out = rho(dv_ind);
    f_new_out = vol_det;
    x_change_out = max(abs(x_new_out-x_old_out));
    f_change_out = max(abs(f_new_out-f_old_out))/(abs(f_new_out)+abs(f_old_out))*2;
    
    fprintf(['OUT Iter.: %d,  Ch.(x): %6.3f,  Ch.(f): %8.5f',repmat('%2d', 1, m_MMA),'\n'], ...
        iter_out,  x_change_out, f_change_out, mpp_iter)
    
    if loopbeta<3, loopbeta = 0; end
end
time_opt = toc(time_start_iter);
fprintf('Finished. Iter.: %3d,  Comp.: %7.3f, Vol.: %7.3f, Ch.: %6.3f, Ch.(f): %8.5f\n', iter, c_mpp, vol_det, x_change_out, f_change_out)


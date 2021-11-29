function pseudograd_rho = pseudograd(par, rho, sx_pseudograd)
% pseudo - gradient computing

pseudograd_rho = zeros(par.dv_num, 2);
for i=1:par.dv_num
    
    num_i = sx_pseudograd.num(i);
    rho_in_sx = rho(sx_pseudograd.nodes_in_sx(1:num_i,i));
    
    rho_without_rhoi = sum(sx_pseudograd.interp(1:num_i, 1) .* rho_in_sx);
    

    pseudograd_rho(i,:) = sum(sx_pseudograd.coeff(1:num_i, i).*(rho_in_sx - rho_without_rhoi).*sx_pseudograd.xj_subs_x(1:num_i,:,i), 1);
end


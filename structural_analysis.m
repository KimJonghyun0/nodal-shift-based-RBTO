function U_xi=structural_analysis(fem, par, sx, rho)

rows = zeros(fem.numelem*36, 1);
cols = zeros(fem.numelem*36, 1);
vals = zeros(fem.numelem*36, 1);
rho_qds = zeros(3,1);
c1 = 1;

for elem = 1:fem.numelem
    for quad_pt = 1:3
        rho_qds(quad_pt) = sx.interp(1:sx.num(elem, quad_pt), elem, quad_pt)'*rho(sx.nodes_in_sx(1:sx.num(elem, quad_pt), elem, quad_pt));
    end
    ke_qd = sum(rho_qds.^par.penal) * fem.bdbwa(:, :, elem);
    edof_tp = fem.edof(elem, :);
    for i1=1:6
        for j1=1:6
            rows(c1) = edof_tp(i1);
            cols(c1) = edof_tp(j1);
            vals(c1) = ke_qd(i1, j1);
            c1=c1+1;
        end
    end
end

K_xi = sparse(rows, cols, vals);
U_xi = zeros(fem.ndof, size(fem.F, 2));

if isfield(fem,'spring')
    K_xi=K_xi+fem.spring;
end
    

U_xi(fem.freedofs, :) = K_xi(fem.freedofs,fem.freedofs)\fem.F(fem.freedofs, :);

% preprocess_findsx
if isempty(cover_coords)
    [~,m,n] = size(xs_quads);
    sx_cover.nodes_in_sx = double.empty(0,m,n);
    sx_cover.dist = double.empty(0,m,n);
    sx_cover.del_xs = double.empty(0,2,m,n);
    sx_cover.num = zeros(m,n);
else
    sx_cover = findsx_direct(xs_quads, cover_coords, par.Rc);
end
sx_cover.nodes_in_sx = sx_cover.nodes_in_sx + par.dv_num;
sx_default = findsx_direct(xs_quads, dv_coords, par.Rc, sx_cover);

tp = findsx_direct(dv_coords, density_coords, par.Rc_pseudograd);
maxtpnum = max(max(tp.num));

nodes_in_sx = ones(maxtpnum-1, par.dv_num);
xj_subs_x = ones(maxtpnum-1, 2, par.dv_num);
coeff = ones(maxtpnum-1, par.dv_num);
interp = ones(maxtpnum-1, par.dv_num);

for i=1:par.dv_num
    num_i = tp.num(i, 1);
    i_pos = find(tp.nodes_in_sx(1:num_i, i)==i);
    nodes_in_sx(1:num_i-1, i) = tp.nodes_in_sx([1:i_pos-1, i_pos+1:num_i], i);
    xj_subs_x(1:num_i-1, :, i) = tp.del_xs([1:i_pos-1, i_pos+1:num_i], :, i);
    dist_sq_tp = tp.dist_sq([1:i_pos-1, i_pos+1:num_i], i);
    f_tp = 1./dist_sq_tp;
    sumf_tp = sum(f_tp);
    coeff(1:num_i-1,i) = -2./(dist_sq_tp.^2) / sumf_tp;
    interp(1:num_i-1,i) = f_tp / sumf_tp;
end

sx_pseudograd.nodes_in_sx = nodes_in_sx;
sx_pseudograd.xj_subs_x = xj_subs_x;
sx_pseudograd.num = tp.num-1;
sx_pseudograd.coeff = coeff;
sx_pseudograd.interp = interp;
ind_sym = zeros(len_cover2, 1);
for i=1:len_cover2
    for j=1:par.dv_num
        if isequal(cover_coords2(i, :).*[1,-1]-dv_coords(j, :), zeros(1, par.dim))
            ind_sym(i) = j;
            break;
        end
    end
end
sx_pseudograd.ind_sym = ind_sym;
clear nodes_in_sx xj_subs_x interp coeff tp dist_tp f_tp df_tp sumf_tp
disp('Sx found');
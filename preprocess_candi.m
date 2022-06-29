% preprocess_candi
max_guess = ceil(4*par.l_candi^2/(par.hx*par.hy));
candi_num = zeros(size(xs_quads, 2), size(xs_quads, 3));
candi_dv = zeros(max_guess, size(xs_quads, 2), size(xs_quads, 3));
for i=1:size(xs_quads,2)
    for j=1:size(xs_quads,3)
        x_diff = xs_quads(:, i, j)' - dv_coords;
        dist_sq = sum(x_diff.^2, 2);

        candi_dv_ij = find(dist_sq<par.l_candi.^2);
        candi_num(i, j) = length(candi_dv_ij);
        candi_dv(1:candi_num(i, j), i, j) = candi_dv_ij;
    end
end
candi_dv = candi_dv(1:max(max(candi_num)), :, :);
disp('Candidates identified');
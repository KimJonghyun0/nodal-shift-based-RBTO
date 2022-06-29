function results = findsx_candi(quad_points, density_points, Rc, dvmax_candi, dxj, candi_dv, candi_num,results_fixed)

[~,m,n] = size(quad_points);
dxj_d = sum(dxj.^2, 2);
more_candi = find(dxj_d > dvmax_candi^2) ;
l_m = length(more_candi);


if nargin == 8
    results_fixed.nodes_in_sx = double.empty(0,m,n);
    results_fixed.del_xs = double.empty(0,2,m,n);
    results_fixed.num = zeros(m,n);    
    results_fixed.dist_sq = double.empty(0,m,n);
    
end

max_guess = 20;
nodes_in_sx = zeros(max_guess, m,n);
del_xs = zeros(max_guess,2, m,n);
dist_sq = zeros(max_guess, m,n);
interp = zeros(max_guess, m,n);
num = zeros(m,n);
denominator = zeros(m,n);
    
    
for i=1:m
    for j=1:n
        % https://stackoverflow.com/questions/50348117/removing-duplicates-in-two-vectors
        i1=1; i2=1; ct=1;
        candi_ij = zeros(candi_num(i, j)+length(more_candi), 1);
        while i1<=candi_num(i, j) && i2<=l_m
            if candi_dv(i1, i, j) < more_candi(i2)
                candi_ij(ct) = candi_dv(i1, i, j);
                ct=ct+1;
                i1=i1+1;
            elseif candi_dv(i1, i, j) > more_candi(i2)
                candi_ij(ct)= more_candi(i2);
                ct=ct+1;
                i2=i2+1;
            else
                candi_ij(ct) = more_candi(i2);
                ct=ct+1;
                i1=i1+1;
                i2=i2+1;
            end
        end
        cc_rem = candi_dv(i1:candi_num(i, j), i, j);
        mc_rem = more_candi(i2:end);
        candi_ij(ct:ct+length(cc_rem)-1) = cc_rem; ct = ct+length(cc_rem);
        candi_ij(ct:ct+length(mc_rem)-1) = mc_rem;
        candi_ij = candi_ij(1:ct+length(mc_rem)-1);
        x_diff = quad_points(:, i, j)' - density_points(candi_ij, :);
        dist_sq_ij = sum(x_diff.^2, 2);
        nodes_in_sx_ij = find(dist_sq_ij - Rc^2 <0);
        num(i,j)  = length(nodes_in_sx_ij) + results_fixed.num(i,j);        
        nodes_in_sx(1:num(i,j), i,j) = [candi_ij(nodes_in_sx_ij); results_fixed.nodes_in_sx(1:results_fixed.num(i,j), i, j)];
        dist_sq(1:num(i,j), i,j) = [dist_sq_ij(nodes_in_sx_ij); results_fixed.dist_sq(1:results_fixed.num(i,j), i,j)];
        del_xs(1:num(i,j), :, i,j) = [x_diff(nodes_in_sx_ij, :); results_fixed.del_xs(1:results_fixed.num(i,j), :, i,j)];
                
        fs = 1./dist_sq(1:num(i,j), i,j);
        denominator(i,j) = sum(fs);
        if denominator(i,j) < inf
            interp(1:num(i, j),i,j) = fs / denominator(i,j);
        else
            tp = zeros(num(i,j), 1); tp(fs==inf) = 1;
            interp(1:num(i, j),i,j) = tp;
        end
    end
end

maxnum = max(max(num));
nodes_in_sx = nodes_in_sx(1:maxnum, :, :);
del_xs = del_xs(1:maxnum, :, :, :);
dist_sq = dist_sq(1:maxnum, :, :);
interp = interp(1:maxnum, :, :);

results.nodes_in_sx = nodes_in_sx;
results.del_xs = del_xs;
results.interp = interp;
results.num = num;
results.dist_sq = dist_sq;
results.denominator = denominator;
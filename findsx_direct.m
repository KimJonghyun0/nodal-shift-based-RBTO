function results = findsx_direct(quad_points, density_points, Rc, intp_func,results_fixed)

[~,m,n] = size(quad_points);
if n==1
    quad_points = quad_points';
    [~,m,n] = size(quad_points);
end

if nargin == 4
    results_fixed.nodes_in_sx = double.empty(0,m,n);
    results_fixed.dist = double.empty(0,m,n);
    results_fixed.del_xs = double.empty(0,2,m,n);
    results_fixed.num = zeros(m,n);    
end


max_guess = 20;
nodes_in_sx = ones(max_guess, m,n);
del_xs = ones(max_guess,2, m,n);
dist = ones(max_guess, m,n);
interp = ones(max_guess, m,n);
num = ones(m,n);
denominator = zeros(m,n);

for i=1:m
    for j=1:n
        x_diff = quad_points(:, i, j)' - density_points;
        
            
        dist_sq = sum(x_diff.^2, 2);
        nodes_in_sx_ij = find(dist_sq - Rc^2 <0);
        num(i,j) = length(nodes_in_sx_ij) + results_fixed.num(i,j);
        nodes_in_sx(1:num(i,j), i,j) = [nodes_in_sx_ij; results_fixed.nodes_in_sx(1:results_fixed.num(i,j), i, j)];
        dist(1:num(i,j), i,j) = [sqrt(dist_sq(nodes_in_sx_ij)); results_fixed.dist(1:results_fixed.num(i,j), i,j)];
        del_xs(1:num(i,j), :, i,j) = [x_diff(nodes_in_sx_ij, :); results_fixed.del_xs(1:results_fixed.num(i,j), :, i,j)];
        
    
        fs = interp_function(dist(1:num(i,j), i,j), intp_func, 0);
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
dist = dist(1:maxnum, :, :);
interp = interp(1:maxnum, :, :);

results.nodes_in_sx = nodes_in_sx;
results.del_xs = del_xs;
results.interp = interp;
results.num = num;
results.dist = dist;
results.denominator = denominator;
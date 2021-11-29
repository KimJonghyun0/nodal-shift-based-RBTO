% Preprocess - FEM preparation


% T3 mesh
try
    switch par.shape
        case 'L-bracket'
            gd = [2, 6, ... % polygon, six line segments
                0, par.lx, par.lx, par.passive_x, par.passive_x, 0, ... % x-coordinates
                0, 0, par.passive_y, par.passive_y, par.ly, par.ly]'; % y-coordinates
            sf = 'P1'; ns = [80;49]; 
        otherwise
            error('Not available yet');
    end
catch % default rectangle
   
end
gd = [3;4;0;par.lx;par.lx;0;0;0;par.ly;par.ly]; ns = char('rect')'; sf='rect';

[dl,bt]=decsg(gd, sf, ns);
model = createpde; pg = geometryFromEdges(model, dl);
generateMesh(model, 'Hmax', 0.5,'GeometricOrder', 'linear');
% pdeplot(model);
fem.nodecoords = model.Mesh.Nodes';
fem.nodenums = model.Mesh.Elements';
fem.numnode = length(fem.nodecoords);
fem.numelem = length(fem.nodenums);
fem.ndof = fem.numnode*par.dim;

% Load
load_vec = []; load_nodes=[];
for i=1:size(load_coords, 1)
    temp = load_coords(i, :);
    if isnan(temp(1))
        load_temp = find(fem.nodecoords(:, 2) == temp(2));
        load_nodes(end+1:end+length(load_temp)) = load_temp;
        load_vec = [load_vec; repmat(load_vector(i, :)/length(load_temp), length(load_temp) , 1)];
    elseif isnan(temp(2))
        load_temp = find(fem.nodecoords(:, 1) == temp(1));
        load_nodes(end+1:end+length(load_temp)) = load_temp;
        load_vec = [load_vec; repmat(load_vector(i, :)/length(load_temp), length(load_temp) , 1)];
        
    else
        load_tp = intersect( find(fem.nodecoords(:, 1) == temp(1)), ...
            find(fem.nodecoords(:, 2) == temp(2)) );
        if isempty(load_tp)
            dist_tp = sum((fem.nodecoords - temp).^2, 2);
            [~, ptn] = min(dist_tp);
            load_nodes(end+1) = ptn;
            warning('load at %.4f, %.4f instead of %.4f, %.4f', ...
                fem.nodecoords(ptn, 1), fem.nodecoords(ptn, 2), ...
                temp(1), temp(2));
        else
            load_nodes(end+1) = intersect( find(fem.nodecoords(:, 1) == temp(1)), ...
                find(fem.nodecoords(:, 2) == temp(2)) );
        end
        load_vec = [load_vec; load_vector(i, :)];
    end
end
fem.F = sparse([load_nodes*2-1, load_nodes*2], 1, load_vec(:), 2*fem.numnode, 1);



% Fixed dofs
fixeddofs = [];
for i=1:size(fixed_coords, 1)
    temp = fixed_coords(i, 1:2);
    if isnan(temp(1))
        fix_temp = find(fem.nodecoords(:, 2) == temp(2));
    elseif isnan(temp(2))
        fix_temp = find(fem.nodecoords(:, 1) == temp(1));
    else
        fix_temp = intersect(find(fem.nodecoords(:, 1) == temp(1)), ...
            find(fem.nodecoords(:, 2) == temp(2)));
    end
    fixeddofs = [fixeddofs; fix_temp*2-2+fixed_coords(i, 3)];
end
alldofs = 1:2*fem.numnode;
fem.freedofs = setdiff(alldofs,fixeddofs);


if exist('compliant_output_coords','var')
    outputdof = [];
    for i=1:size(compliant_output_coords, 1)
        temp = compliant_output_coords(i, 1:2);
        if isnan(temp(1))
            comp_temp = find(fem.nodecoords(:, 2) == temp(2));
        elseif isnan(temp(2))
            comp_temp = find(fem.nodecoords(:, 1) == temp(1));
        else
            comp_temp = intersect(find(fem.nodecoords(:, 1) == temp(1)), ...
                find(fem.nodecoords(:, 2) == temp(2)));
        end
        outputdof = [outputdof; comp_temp*2-2+compliant_output_coords(i, 3)];
    end
    fem.F(outputdof, 2) = 1;
    fem.outputdof = outputdof;
end

if exist('spring_coords', 'var')
    springdof = []; springcoeff = [];
    for i=1:size(spring_coords, 1)
        temp = spring_coords(i, 1:2);
        if isnan(temp(1))
            comp_temp = find(fem.nodecoords(:, 2) == temp(2));
        elseif isnan(temp(2))
            comp_temp = find(fem.nodecoords(:, 1) == temp(1));
        else
            comp_temp = intersect(find(fem.nodecoords(:, 1) == temp(1)), ...
                find(fem.nodecoords(:, 2) == temp(2)));
        end
        springdof = [springdof; comp_temp*2-2+spring_coords(i, 3)];
        springcoeff = [springcoeff; spring_coords(i, 4) * ones(size(comp_temp, 1), 1)];
    end
    fem.spring = sparse(springdof, springdof, springcoeff, fem.ndof, fem.ndof);
end


fem.D0 = D0;
fem.weights_qd=1/3;




% Quadrature points
area = -ones(fem.numelem, 1);
xs_quads = ones(2, fem.numelem,3);
bdbwa = ones(6, 6, fem.numelem);
edof = -ones(fem.numelem, par.dim*3);
xs_quads_matrix = zeros(fem.numelem*3, 2);
for elem = 1:fem.numelem
    x_t3 = fem.nodecoords(fem.nodenums(elem, :), 1);
    y_t3 = fem.nodecoords(fem.nodenums(elem, :), 2);
    quadcoords(:, 1) = ones(3, 1) * sum(x_t3)/6;
    quadcoords(:, 2) = ones(3, 1) * sum(y_t3)/6;
    quadcoords = quadcoords + [x_t3, y_t3]/2;
    
    for quad_pt = 1:3
        x_qdpt = quadcoords(quad_pt, :);  %
        xs_quads(1:2, elem, quad_pt) = x_qdpt;
        xs_quads_matrix(elem*3-3+quad_pt, :) = x_qdpt;
    end
    area(elem) = abs(det([x_t3';y_t3';[1,1,1]])/2);
    
    
    den_Bqd = x_t3(2)*y_t3(3) - x_t3(3)*y_t3(2) + x_t3(3)*y_t3(1) - x_t3(1)*y_t3(3) + x_t3(1)*y_t3(2) - x_t3(2)*y_t3(1);
    B_qd = [y_t3(2)-y_t3(3)        0        y_t3(3)-y_t3(1)        0        y_t3(1)-y_t3(2)        0;
        0           x_t3(3)-x_t3(2)        0        x_t3(1)-x_t3(3)        0        x_t3(2)-x_t3(1);
        x_t3(3)-x_t3(2) y_t3(2)-y_t3(3) x_t3(1)-x_t3(3) y_t3(3)-y_t3(1) x_t3(2)-x_t3(1) y_t3(1)-y_t3(2)]/den_Bqd;
    bdbwa(:, :, elem) = B_qd'*D0*B_qd*fem.weights_qd*area(elem); % B_qd'*D0*B_qd*weight_qd
    
    edof(elem, :) = [2*fem.nodenums(elem, 1) + [-1, 0], 2*fem.nodenums(elem, 2) ...
        + [-1, 0], 2*fem.nodenums(elem, 3) + [-1, 0]];
end

fem.edof = edof;
fem.bdbwa = bdbwa;
fem.area = area;

clear edof bdbwa area

disp('Domain meshed');

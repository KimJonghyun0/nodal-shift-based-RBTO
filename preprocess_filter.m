if ~exist('r_filter', 'var'), r_filter = 2; end
% PREPARE FILTER
iH = zeros((par.dv_num+len_cover)*(par.dv_num+len_cover),1); iH=cast(iH, 'uint64');
jH = zeros(size(iH)); jH=cast(jH, 'uint64');
sH = zeros(size(iH));
k = 0;         
for i=1:par.dv_num+len_cover
    for j=1:par.dv_num+len_cover
        k = k+1;
        iH(k) = i;
        jH(k) = j;
        sH(k) = max(0, r_filter - sqrt( sum((density_coords(i, :)-density_coords(j, :)).^2)));
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
clearvars iH jH
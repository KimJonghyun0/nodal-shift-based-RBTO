% preprocess - Karhunen Loeve expansion
if dtomode == 1
    kl = 0;
    disp('Passing K-L expansion: Unnecessary for DTO')
    return;
end

kl.phi = []; kl.mean = [];
for i_dim = 1:size(sigma_corr, 2)
    r=1:(par.nely+1)*(par.nelx+1); c=r';
    cov=sigma_corr(r, i_dim)'.*sigma_corr(c, i_dim).*exp( ...
        -((par.hx*(floor((r-1)/(par.nely+1))-floor((c-1)/(par.nely+1)))).^2 ...
        +(par.hy*(r-(par.nely+1)*(floor((r-1)/(par.nely+1)))-c+(par.nely+1)*(floor((c-1)/(par.nely+1))))).^2) ...
        /(par.corr_length^2));   
    b=par.hx/2; h = par.hy/2;
    [eigvec, eigval] = eigs(cov, par.eignum);
    eigval = eigval * par.hx * par.hy;
    err = 1;    
    intg_sigsq = sum(sigma_corr(:, i_dim).^2).*par.hx*par.hy;   
    i=0; phi = [];
    while err > error_target && i < par.eignum
        i=i+1;
        err = err - eigval(i,i)/intg_sigsq;
        
        phi_i = eigvec(:, i)' * cov/sqrt(eigval(i,i));
        
        phi = [phi, phi_i(:)];
    end
    % error
    if err > error_target
        warning('Need more eigenvalues');
        par.eignum_extended = par.eignum;
        while err>error_target
            i = par.eignum_extended;
            if par.eignum_extended >= length(cov)
                err('par.eignum_extended >= node numbers in K-L phase')
            end
            par.eignum_extended = min(2*par.eignum_extended, length(cov));
            [eigvec, eigval] = eigs(cov, par.eignum_extended);
            while err > error_target && i < par.eignum_extended
                i=i+1;
                err = err - eigval(i,i)/intg_sigsq;
                phi_i = eigvec(:, i)' * cov/sqrt(eigval(i,i));
                phi = [phi, phi_i(:)];
            end
        end
    end
    kl.m(i_dim) = i;
    kl.phi = [kl.phi, phi(:, 1:i)];
    try
        kl.mean = [kl.mean, mean_corr];
    catch
        kl.mean = [kl.mean, zeros(par.dv_num, 1)];
    end
    clear cov    
end
kl.cumsum_m = [0, cumsum(kl.m)];
kl.dim = size(sigma_corr, 2);
disp('K-L expansion done');
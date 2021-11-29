% preprocess - Karhunen Loeve expansion
if dtomode == 1
    kl = 0; 
    disp('Passing K-L expansion: Unnecessary for DTO')
    return; 
end


r=1:(par.nely+1)*(par.nelx+1); c=r';
cov=sigma_corr(r)'.*sigma_corr(c).*exp( ...
    -((par.hx*(floor((r-1)/(par.nely+1))-floor((c-1)/(par.nely+1)))).^2 ...
+(par.hy*(r-(par.nely+1)*(floor((r-1)/(par.nely+1)))-c+(par.nely+1)*(floor((c-1)/(par.nely+1))))).^2) ...
    /(par.corr_length^2));


b=par.hx/2; h = par.hy/2;
[eigvec, eigval] = eigs(cov, par.eignum); 
eigval = eigval * par.hx * par.hy; 
err = 1;

intg_sigsq=0;
intg_ntn = [(4*b*h)/9, (2*b*h)/9,   (b*h)/9, (2*b*h)/9;
 (2*b*h)/9, (4*b*h)/9, (2*b*h)/9,   (b*h)/9;
   (b*h)/9, (2*b*h)/9, (4*b*h)/9, (2*b*h)/9;
 (2*b*h)/9,   (b*h)/9, (2*b*h)/9, (4*b*h)/9];
for ely = 1:nely
    for elx = 1:nelx
        n1 = (nely+1)*(elx-1)+ely;
        n2 = (nely+1)* elx   +ely;
        esigma = sigma_corr([n1; n2; n2+1;n1+1]);
        intg_sigsq = intg_sigsq + esigma' * intg_ntn * esigma;
    end
end

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
        i = par.eignum_extended + 1;
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


kl.m = i;
kl.phi = phi(:, 1:kl.m);
kl.ev = diag(eigval(1:kl.m, 1:kl.m));
try 
    kl.mean = mean_corr;
catch
    kl.mean = zeros(par.dv_num, 1); 
end
clear cov

disp('K-L expansion done');
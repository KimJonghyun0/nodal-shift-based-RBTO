function [m_MMA, n_MMA, iter, xval, xold1, xold2, low, upp, xmin, xmax, a0_MMA, a_MMA, c_MMA, d_MMA, MMA_movelimit] =...
    initialize_MMA(dv_num, rhomin)
m_MMA = 1; n_MMA = dv_num; iter = 0; 
xval = 0.5*ones(n_MMA, 1);
xold1 =xval; xold2= xval; low = xval; upp = xval;
xmin=rhomin*ones(n_MMA, 1); xmax=ones(n_MMA, 1);
a0_MMA = 1; a_MMA = zeros(m_MMA,1); c_MMA = 1.0*10^6*ones(m_MMA,1); d_MMA = ones(m_MMA,1);
MMA_movelimit = 1-rhomin;
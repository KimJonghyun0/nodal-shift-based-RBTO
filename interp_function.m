function fdf=interp_function(d, intp_func, order)
if intp_func == 1 % Shepard
    if order == 0
        fdf = 1./d.^2;
    elseif order == 1
        fdf = -2./d.^3;
    end
elseif intp_func == 2 % exp(-d)
    if order == 0
        fdf = exp(-d);
    elseif order == 1
        fdf = -exp(-d);
    end
elseif intp_func == 3 % exp(-d^2)
    if order == 0
        fdf = exp(-d.^2);
    elseif order == 1
        fdf = -2*d.*exp(-d.^2);
    end
elseif intp_func == 4 % Shepard_1
    if order == 0
        fdf = 1./d;
    elseif order == 1
        fdf = -1./d.^2;
    end
elseif intp_func == 5 % Shepard_3
    if order == 0
        fdf = 1./d.^3;
    elseif order == 1
        fdf = -3./d.^4;
    end
elseif intp_func == 6 % Shepard_5
    if order == 0
        fdf = 1./d.^5;
    elseif order == 1
        fdf = -5./d.^6;
    end
elseif intp_func == 7 % exp(-d^3)
    if order == 0
        fdf = exp(-d.^3);
    elseif order == 1
        fdf = -3*d.^2.*exp(-d.^3);
    end
elseif intp_func == 8 % exp(1/d)-1
    if order == 0
        fdf = exp(1./d)-1;
    elseif order == 1
        fdf = -1./d.^2.*exp(1./d);
    end
elseif intp_func == 9 % exp(1/d^2)-1
    if order == 0
        fdf = exp(1./d.^2)-1;
    elseif order == 1
        fdf = -2./d.^3.*exp(1./d.^2);
    end
elseif intp_func == 10 % exp(1/d^3)-1
    if order == 0
        fdf = exp(1./d.^3)-1;
    elseif order == 1
        fdf = -3./d.^4.*exp(1./d.^3);
    end
elseif intp_func == 11 % 1/log(d+1)
    if order == 0
        fdf = 1./log(d+1);
    elseif order == 1
        fdf = -1./((log(d+1)).^2 .* (d+1));
    end
elseif intp_func == 12 % 1/log(d+1)^2
    if order == 0
        fdf = 1./log(d+1).^2;
    elseif order == 1
        fdf = -2./((log(d+1)).^3 .* (d+1));
    end
elseif intp_func == 13 % 1/log(d+1)^3
    if order == 0
        fdf = 1./log(d+1).^3;
    elseif order == 1
        fdf = -3./((log(d+1)).^4 .* (d+1));
    end
elseif intp_func == 14 %  exp(1/d^0.2)-1
    if order == 0
        fdf = exp(1./d.^0.2)-1;
    elseif order == 1
        fdf = -0.2./d.^1.2.*exp(1./d.^0.2);
    end
elseif intp_func == 15 %  exp(1/d^0.1)-1
    if order == 0
        fdf = exp(1./d.^0.1)-1;
    elseif order == 1
        fdf = -0.1./d.^1.1.*exp(1./d.^0.1);
    end
elseif intp_func == 16 %  csch
    if order == 0
        fdf = csch(d);
    elseif order == 1
        fdf = -coth(d).*csch(d);
    end
elseif intp_func == 17 %  coth
    if order == 0
        fdf = coth(d) - 1;
    elseif order == 1
        fdf = -(csch(d)).^2;
    end
else
    error('Undefined interpolation function')
end
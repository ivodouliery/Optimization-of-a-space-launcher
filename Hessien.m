function [hess] = Hessien(hess_m1, gf, gc, gf_m1, gc_m1, xk, xk_m1, lambdak, init)
n = length(xk_m1);
if init == 0
    hess = eye(n);
else
    gL_m1 = gf_m1 + gc_m1'*lambdak;
    gL = gf + gc'*lambdak;

    dk_m1 = xk - xk_m1;
    yk_m1 = gL - gL_m1;

    if init == 1
        if yk_m1' * dk_m1 > 0
            hess = hess_m1 + (yk_m1*yk_m1')/(yk_m1'*dk_m1) - (hess_m1*dk_m1*(dk_m1')*hess_m1)/(dk_m1'*hess_m1*dk_m1);
        else
            hess = hess_m1;
        end
        
    elseif init == 2
        if dk_m1'*(yk_m1-hess_m1*dk_m1) ~= 0
            hess = hess_m1 + ((yk_m1-hess_m1*dk_m1)*(yk_m1-hess_m1*dk_m1)')/(dk_m1'*(yk_m1-hess_m1*dk_m1));
        else
            hess = hess_m1;
        end
    end
end
end
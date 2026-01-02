function [xk_p1, nb_f, merite_kp1, n_fonc] = Globalisation(fc, xk, d_QP, rho, nb_fail, nfonc, cx, gf)

d_merite = gf'*d_QP - rho * norm(cx);
[merite_k, n_fonc] = Merite(fc, xk, rho, nfonc);

if d_merite < 1e-10
    s=1;
    xk_p1 = xk+s*d_QP;
    [merite_kp1, n_fonc] = Merite(fc, xk_p1, rho, n_fonc);
    cpt = 0;

    while merite_kp1 >= merite_k + 0.1 * s * d_merite
        s = s/2;
        xk_p1 = xk+s*d_QP;
        [merite_kp1, n_fonc] = Merite(fc, xk_p1, rho, n_fonc);
        if cpt == 10
            break;
        end
        cpt = cpt+1;
    end

    nb_f = 0;
else
    xk_p1 = xk;
    nb_f = nb_fail + 1;
    merite_kp1 = merite_k;
end
if nb_f > 5
    error('Globalisation échouée après 5 tentatives');
end
end
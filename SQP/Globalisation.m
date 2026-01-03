function [xk_p1, nb_f, merite_kp1, n_fonc, fx_p1, cx_p1] = Globalisation(fc, xk, d_QP, rho, nb_fail, nfonc, cx, fx, gf, c1, lb, ub)
% Globalisation Recherche lineaire (Armijo)

d_merite = gf'*d_QP - rho * norm(cx, 1);
merite_k = fx + rho*norm(cx, 1);

if d_merite < 1e-10
    s=1;
    xk_p1 = xk+s*d_QP;

    % Clamping (Bornes)
    xk_p1 = max(lb, min(ub, xk_p1));

    % Init
    [fx_p1, cx_p1] = fc(xk_p1);
    n_fonc = nfonc+1;
    merite_kp1 = fx_p1 + rho*norm(cx_p1, 1);

    cpt = 0;
    while merite_kp1 >= merite_k + c1 * s * d_merite
        s = s/2;
        xk_p1 = xk+s*d_QP;

        % Clamping updates
        xk_p1 = max(lb, min(ub, xk_p1));

        % Eval
        [fx_p1, cx_p1] = fc(xk_p1);
        n_fonc = nfonc+1;
        merite_kp1 = fx_p1 + rho*norm(cx_p1, 1);

        if cpt == 20
            break;
        end
        cpt = cpt+1;
    end
    nb_f = 0;
else
    xk_p1 = xk;
    fx_p1 = fx;
    cx_p1 = cx;
    nb_f = nb_fail + 1;
    merite_kp1 = merite_k;
    n_fonc = nfonc;
end

if nb_f > 20
    error('Echec Globalisation');
end
end
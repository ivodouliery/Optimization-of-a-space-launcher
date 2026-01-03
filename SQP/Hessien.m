function [hess] = Hessien(hess_m1, gf, gc, gf_m1, gc_m1, xk, xk_m1, lambdak, init)
% Mise à jour du Hessien (BFGS ou SR1)

n = length(xk_m1);
if init == 0
    hess = eye(n); % Identité au départ
else
    % Gradient du Lagrangien
    gL_m1 = gf_m1 + gc_m1'*lambdak;
    gL = gf + gc'*lambdak;

    dk_m1 = xk - xk_m1;   % Déplacement x
    yk_m1 = gL - gL_m1;   % Variation gradient

    if init == 1
        % Mise à jour BFGS
        if yk_m1' * dk_m1 > 1e-8
            hess = hess_m1 + (yk_m1*yk_m1')/(yk_m1'*dk_m1) - (hess_m1*dk_m1*(dk_m1')*hess_m1)/(dk_m1'*hess_m1*dk_m1);
        else
            hess = hess_m1;
        end

    elseif init == 2
        % Mise à jour SR1
        if abs(dk_m1'*(yk_m1-hess_m1*dk_m1)) > 1e-8
            hess = hess_m1 + ((yk_m1-hess_m1*dk_m1)*(yk_m1-hess_m1*dk_m1)')/(dk_m1'*(yk_m1-hess_m1*dk_m1));
        else
            hess = hess_m1;
        end
    end
end
end
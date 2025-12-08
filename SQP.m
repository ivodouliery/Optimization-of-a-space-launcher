function [x, lambda] = SQP(fc, x0, lambda0, h, mod)

% Point initial
xk = x0;
lambdak = lambda0;

k = 0;
init = 0;
nb_fail = 0;
Q = 0;

[~, ck] = fc(xk);
n = length(x0);
m = length(ck);

gf_m1 = zeros(n, 1);
gc_m1 = zeros(m, n);
xk_m1 = zeros(n, 1);
merite_m1 = 0;

rho = 1;

k_max = 50;
tol_x = 1e-10;
tol_f = 1e-10;
nfonc = 1;

while k<k_max % Tests d'arrêt
    % Gradient
    [gf, gc, nfonc] = gradient(fc,xk,h, nfonc);

    % Quasi-Newton
    if nb_fail > 0
        init = 0;
    end
    Q = Hessien(Q, gf, gc, gf_m1, gc_m1, xk, xk_m1, lambdak, init);
    init = mod;

    % Modification de H
    Q = modif_hessien(Q);

    % Solution problème quadratique
    [lambdak, d_QP] = Sol_pb_quadra(Q, gf, gc, ck);

    % Globalisation
    if nb_fail > 1 
        rho = rho * 10;
    end
    [xk, nb_fail, merite_k, nfonc] = Globalisation(fc, xk, d_QP, rho, nb_fail, nfonc);
    
    %Test d'arrêt sur déplacement de x
    if norm(xk-xk_m1) < tol_x && k > 0 && nb_fail == 0
        gradL = gf + sum(lambdak .* gc);
        affiche_iteration(k, nfonc, xk, fk, ck, lambdak, gradL);
        disp('Arrêt : Déplacement de x insuffisant');
        break; 
    end

    [fk, ck] = fc(xk);
    nfonc = nfonc+1;

    if nb_fail == 0

        %Test d'arrêt sur amélioration de f
        if abs(merite_k - merite_m1) < tol_f && k > 0
            gradL = gf + gc'*lambdak;
            affiche_iteration(k, nfonc, xk, fk, ck, lambdak, gradL);
            disp('Arrêt : Déplacement de f insuffisant');
            break; 
        end

        [gf_m1, gc_m1] = deal(gf, gc);
        xk_m1 = xk;
        merite_m1 = merite_k;

        gradL = gf + sum(lambdak .* gc);
        affiche_iteration(k, nfonc, xk, fk, ck, lambdak, gradL);

        k = k+1;
    end

    
    if k == k_max
        disp('Arrêt : itérations max atteint')
    end

end



x = xk;
lambda = lambdak;
end
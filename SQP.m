function [x, lambda] = SQP(fc, x0, lambda0, lb, ub, options)

% Extraction des options
if isfield(options, 'step_size'), h = options.step_size; else, h = 1e-7; end
if isfield(options, 'hessian_mod'), mod = options.hessian_mod; else, mod = 1; end
if isfield(options, 'max_iter'), k_max = options.max_iter; else, k_max = 50; end
if isfield(options, 'max_eval'), max_eval = options.max_eval; else, max_eval = 1000; end
if isfield(options, 'tol_x'), tol_x = options.tol_x; else, tol_x = 1e-10; end
if isfield(options, 'tol_f'), tol_f = options.tol_f; else, tol_f = 1e-10; end

% Point initial
xk = x0;
lambdak = lambda0;

k = 0;
init = 0;
nb_fail = 0;
Q = 0;

[fk, ck] = fc(xk);
n = length(x0);
m = length(ck);

gf_m1 = zeros(n, 1); % Permet de stocker le gradient précédent
gc_m1 = zeros(m, n); % Permet de stocker le gradient des contraintes précédent
xk_m1 = zeros(n, 1); % Permet de stocker le point précédent
merite_m1 = 0; % Permet de stocker le merite précédent

rho = 1;
nfonc = 1;

while k < k_max && nfonc < max_eval % Tests d'arrêt
    % Gradient
    [gf, gc, nfonc] = gradient(fc,xk,h, nfonc, fk, ck);

    % Quasi-Newton
    if nb_fail > 0
        init = 0;
    end
    Q = Hessien(Q, gf, gc, gf_m1, gc_m1, xk, xk_m1, lambdak, init);
    init = mod;

    % Modification de H
    Q = modif_hessien(Q);
    a
    % Solution problème quadratique
    [lambdak, d_QP] = Sol_pb_quadra(Q, gf, gc, ck);

    % Globalisation
    if nb_fail > 1
        rho = rho * 10;
    end
    [xk, nb_fail, merite_k, nfonc] = Globalisation(fc, xk, d_QP, rho, nb_fail, nfonc, ck, gf);

    %Test d'arrêt sur déplacement de x
    if norm(xk-xk_m1) < tol_x && k > 0 && nb_fail == 0
        gradL = gf + gc' * lambdak;
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

        gradL = gf + gc' * lambdak;
        affiche_iteration(k, nfonc, xk, fk, ck, lambdak, gradL);

        k = k+1;
    end


    if k == k_max
        disp('Arrêt : itérations max atteint')
    elseif nfonc >= max_eval
        disp('Arrêt : nombre max d évaluations atteint')
    end

end



x = max(lb, min(ub, xk));
lambda = lambdak;
end
function [x, lambda] = SQP(fc, x0, lambda0, lb, ub, options)
% Algorithme SQP (Sequential Quadratic Programming) avec bornes

fprintf('=== SQP ===\n');

% --- Options et Initialisation ---
if isfield(options, 'step_size'), h = options.step_size; else, h = 1e-6 * max(abs(x0), 1.0); end
if isfield(options, 'hessian_mod'), mod = options.hessian_mod; else, mod = 1; end
if isfield(options, 'max_iter'), k_max = options.max_iter; else, k_max = 100; end
if isfield(options, 'max_eval'), max_eval = options.max_eval; else, max_eval = 1000; end
if isfield(options, 'tol_x'), tol_x = options.tol_x; else, tol_x = 1e-10; end
if isfield(options, 'tol_f'), tol_f = options.tol_f; else, tol_f = 1e-10; end

xk = x0;
lambdak = lambda0;
k = 0;
nb_fail = 0;
Q = 0; % Hessien initial (sera Id)
nfonc = 1;
rho = 1; % Paramètre de pénalité

[fk, ck] = fc(xk);
n = length(x0);
m = length(ck);

% Variables pour mémoire de l'étape précédente
gf_m1 = zeros(n, 1);
gc_m1 = zeros(m, n);
xk_m1 = zeros(n, 1);
merite_m1 = 0;

% --- Boucle Principale ---
while k < k_max && nfonc < max_eval
    % 1. Calcul des Gradients
    [gf, gc, nfonc] = gradient(fc,xk,h, nfonc, fk, ck);

    % 2. Mise à jour du Hessien (BFGS/SR1)
    if nb_fail > 0, init = 0; else, init = mod; end % Reset si echec
    if k==0, init=0; end
    Q = Hessien(Q, gf, gc, gf_m1, gc_m1, xk, xk_m1, lambdak, init);

    % 3. Régularisation (si non défini positif)
    Q = modif_hessien(Q);

    % 4. Résolution du Sous-Problème Quadratique (QP) -> Direction d
    [lambdak, d_QP] = Sol_pb_quadra(Q, gf, gc, ck);

    % 5. Gestion paramètre de pénalité (Augmentation si difficultés)
    if nb_fail > 1, rho = rho * 10; end

    % Critère Armijo strict ou relâché
    if k < k_max/2, c1 = 1e-3; else, c1 = 1e-1; end

    % 6. Recherche Linéaire (Globalisation)
    [xk, nb_fail, merite_k, nfonc, fk_new, ck_new] = Globalisation(fc, xk, d_QP, rho, nb_fail, nfonc, ck, fk, gf, c1, lb, ub);

    % 7. Projection finale sur bornes (Sécurité)
    xk = max(lb, min(ub, xk));

    % --- Tests d'arrêt ---
    % Sur déplacement x
    if norm(xk - xk_m1) < tol_x && k > 0 && nb_fail == 0
        affiche_iteration(k, nfonc, xk, fk, ck, lambdak, gf + gc' * lambdak);
        disp('Arrêt : Convergence x');
        break;
    end

    [fk, ck] = deal(fk_new, ck_new);

    if nb_fail == 0
        % Sur amélioration fonction mérite
        if abs(merite_k - merite_m1) < tol_f && k > 0
            affiche_iteration(k, nfonc, xk, fk, ck, lambdak, gf + gc'*lambdak);
            disp('Arrêt : Stagnation');
            break;
        end

        % Sauvegarde pour itération suivante
        [gf_m1, gc_m1] = deal(gf, gc);
        xk_m1 = xk;
        merite_m1 = merite_k;

        % Affichage
        gradL = gf + gc' * lambdak;
        affiche_iteration(k, nfonc, xk, fk, ck, lambdak, gradL);

        k = k+1;
    end

    if k == k_max, disp('Arrêt : Max itérations');
    elseif nfonc >= max_eval, disp('Arrêt : Max évaluations'); end

end

x = max(lb, min(ub, xk));
lambda = lambdak;
end
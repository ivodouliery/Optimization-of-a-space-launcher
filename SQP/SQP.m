function [x, lambda] = SQP(fc, x0, lambda0, lb, ub, options)
% SQP Algorithme de programmation quadratique sequentielle

fprintf('=== SQP ===\n');

% --- Options ---
if isfield(options, 'step_size'), h = options.step_size; else, h = 1e-6 * max(abs(x0), 1.0); end
if isfield(options, 'hessian_mod'), mod = options.hessian_mod; else, mod = 1; end
if isfield(options, 'max_iter'), k_max = options.max_iter; else, k_max = 100; end
if isfield(options, 'max_eval'), max_eval = options.max_eval; else, max_eval = 1000; end
if isfield(options, 'tol_x'), tol_x = options.tol_x; else, tol_x = 1e-10; end
if isfield(options, 'tol_f'), tol_f = options.tol_f; else, tol_f = 1e-10; end

% --- Init ---
xk = x0;
lambdak = lambda0;
k = 0;
init = 0;
nb_fail = 0;
Q = 0;

[fk, ck] = fc(xk);
n = length(x0);
m = length(ck);

% Memoire
gf_m1 = zeros(n, 1);
gc_m1 = zeros(m, n);
xk_m1 = zeros(n, 1);
merite_m1 = 0;

rho = 1;
nfonc = 1;

% --- Boucle ---
while k < k_max && nfonc < max_eval
    % Gradient
    [gf, gc, nfonc] = gradient(fc,xk,h, nfonc, fk, ck);

    % Hessien
    if nb_fail > 0, init = 0; end
    Q = Hessien(Q, gf, gc, gf_m1, gc_m1, xk, xk_m1, lambdak, init);
    init = mod;
    Q = modif_hessien(Q);

    % QP
    [lambdak, d_QP] = Sol_pb_quadra(Q, gf, gc, ck);

    % Globalisation
    if nb_fail > 1, rho = rho * 10; end
    if k < k_max/2, c1 = 1e-3; else, c1 = 1e-1; end

    [xk, nb_fail, merite_k, nfonc, fk_new, ck_new] = Globalisation(fc, xk, d_QP, rho, nb_fail, nfonc, ck, fk, gf, c1, lb, ub);

    % Bornes
    xk = max(lb, min(ub, xk));

    % Arret (x)
    if norm(xk - xk_m1) < tol_x && k > 0 && nb_fail == 0
        gradL = gf + gc' * lambdak;
        affiche_iteration(k, nfonc, xk, fk, ck, lambdak, gradL);
        disp('Arret : dx < tol');
        break;
    end

    [fk, ck] = deal(fk_new, ck_new);

    if nb_fail == 0
        % Arret (f)
        if abs(merite_k - merite_m1) < tol_f && k > 0
            gradL = gf + gc'*lambdak;
            affiche_iteration(k, nfonc, xk, fk, ck, lambdak, gradL);
            disp('Arret : df < tol');
            break;
        end

        % Update memoire
        [gf_m1, gc_m1] = deal(gf, gc);
        xk_m1 = xk;
        merite_m1 = merite_k;

        gradL = gf + gc' * lambdak;
        affiche_iteration(k, nfonc, xk, fk, ck, lambdak, gradL);
        k = k+1;
    end

    if k == k_max
        disp('Arret : Max iter');
    elseif nfonc >= max_eval
        disp('Arret : Max eval');
    end
end

x = max(lb, min(ub, xk));
lambda = lambdak;
end
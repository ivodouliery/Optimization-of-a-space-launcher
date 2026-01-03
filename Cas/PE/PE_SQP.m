clear;
close all;

k = [0.1; 0.15; 0.2];
Ve = [2600; 3000; 4400];
Vp = 10000;
mu_sat = 1000;

x0 = [3; 3; 3];
lambda0 = 0;
lb = [0.1; 0.1; 0.1];
ub = [10.9; 7.6; 5.9];

options.hessian_mod = 2;
options.max_iter = 1000;
options.max_eval = 5000;
options.tol_x = 1e-8;
options.tol_f = 1e-8;

fprintf('Lancement optimisation PE (Ratios x_i)...\n');

% Definition de la fonction de scaling locale
fc = @(x) f_PE(x, k, Ve, Vp);

[x_star, lambda_star] = SQP(fc, x0, lambda0, lb, ub, options);

fprintf('\nRatios optimaux trouvés :\n');
fprintf('x1 = %.4f\n', x_star(1));
fprintf('x2 = %.4f\n', x_star(2));
fprintf('x3 = %.4f\n', x_star(3));

% Calcul des masses physiques pour validation
Mi_next = mu_sat;
masses_ergols = zeros(3,1);
for j = 3:-1:1
    y_j = (1 + k(j)) / x_star(j) - k(j);
    Mi_curr = Mi_next / y_j;

    % Masse ergol : me = Mi * (x-1)/x
    masses_ergols(j) = Mi_curr * (x_star(j) - 1) / x_star(j);

    Mi_next = Mi_curr; % Pour l'étage du dessous
end
M0_calc = Mi_next;

fprintf('\nMasses d''ergols calculées (kg) :\n');
fprintf('me1 = %.4f kg\n', masses_ergols(1));
fprintf('me2 = %.4f kg\n', masses_ergols(2));
fprintf('me3 = %.4f kg\n', masses_ergols(3));

[f_val, ~] = fc(x_star);
ratio_J = -f_val;

fprintf('\nPerformance :\n');
fprintf('Ratio Charge Utile J = %.6f\n', ratio_J);
fprintf('Masse totale M0 (via J) = %.4f kg\n', mu_sat / ratio_J);



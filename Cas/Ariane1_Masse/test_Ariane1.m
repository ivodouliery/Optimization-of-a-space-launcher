clear;
close all;

scale = 10000;
x0 = [150000; 40000; 10000] / scale;
lambda0 = 0;
lb = [1000; 1000; 100] / scale;
ub = [1e6; 1e6; 1e6] / scale;

options.step_size = 1e-4; % Adjusted for scaled variables (order 10)
options.hessian_mod = 1;
options.max_iter = 100;
options.max_eval = 2000;
options.tol_x = 1e-6;
options.tol_f = 1e-6;

fprintf('Lancement optimisation Ariane 1\n');

% Definition de la fonction de scaling locale
f_scaled = @(x) f_wrapper(x, scale);

[x_star_scaled, lambda_star] = SQP(f_scaled, x0, lambda0, lb, ub, options);

x_star = x_star_scaled * scale;

fprintf('\nRésultat trouvé :\n');
fprintf('me1 = %.4f kg\n', x_star(1));
fprintf('me2 = %.4f kg\n', x_star(2));
fprintf('me3 = %.4f kg\n', x_star(3));

[M0, ~] = f_Ariane1(x_star);
fprintf('Masse totale M0 = %.4f kg\n', M0);


function [f, c] = f_wrapper(x, scale)
x_real = x * scale;
[f_real, c] = f_Ariane1(x_real);

% Scale objective to be order 1 (approx)
f = f_real / scale;
end

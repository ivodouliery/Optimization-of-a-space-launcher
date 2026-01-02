clear;
close all;


x0 = [-1;2;1;-2;-2];
lambda0 = zeros(3,1);

ub = [0; 3; 2; 0; -1];
lb = [-2; 1; 0; -3; -3];

options.step_size = 1e-10;
options.hessian_mod = 1;
options.max_iter = 50;
options.max_eval = 1000;
options.tol_x = 1e-4;
options.tol_f = 1e-4;

[x_star, lambda_star] = SQP(@f_MHW4D, x0, lambda0, lb, ub, options);
fprintf('\nRésultat trouvé :\n');
fprintf('x_1 = %1.4f\n', x_star(1));
fprintf('x_2 = %1.4f\n', x_star(2));
fprintf('x_3 = %1.4f\n', x_star(3));
fprintf('x_4 = %1.4f\n', x_star(4));
fprintf('x_5 = %1.4f\n', x_star(5));
[f_star, ~] = f_MHW4D(x_star);
fprintf('f_star = %1.4f\n', f_star);




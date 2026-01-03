function [] = affiche_iteration(k, nfonc, xk, fxk, cxk, lambdak, gradL_k)

fprintf('iter: %-4d nfonc: %-6d f(x): %-12.4f ||c(x)||: %-12.4f ||gradL||: %-12.4f\n', ...
    k, nfonc, fxk, norm(cxk), norm(gradL_k));
fprintf('    xk:     %s\n', sprintf('%12.4f', xk));
fprintf('    lambdak: %s\n\n', sprintf('%12.4f', lambdak));

end
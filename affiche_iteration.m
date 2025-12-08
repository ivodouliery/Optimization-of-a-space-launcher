function [] = affiche_iteration(k, nfonc, xk, fxk, cxk, lambdak, gradL_k)


    xk_str = sprintf('%.4f, ', xk);
    xk_str = ['(' xk_str(1:end-2) ')'];

    lam_str = sprintf('%.4f, ', lambdak);
    lam_str = ['(' lam_str(1:end-2) ')'];

    fprintf('iter: %d\t nfonc: %d\t xk: %s\t f(x): %.4f\t c(x): %.4f\t lambda: %s\t ||gradL||: %.4f\n', k, nfonc, xk_str, fxk, norm(cxk), lam_str, norm(gradL_k));
end

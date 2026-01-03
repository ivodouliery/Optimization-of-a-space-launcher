function [gf, gc, n_fonc] = gradient(fc, x, h, nfonc, fx, cx)
% Calcul du gradient par différences finies

fk = fx;
ck = cx;
n_fonc = nfonc;
n = length(x);
m = length(ck);
gf = zeros(n,1); % Gradient de la fonction objectif f
gc = zeros(m,n); % Jacobienne des contraintes c

for i=1:n
    xk = x;
    % Pas de perturbation adapté
    if isscalar(h), hi = h; else, hi = h(i); end
    xk(i) = x(i) + hi;

    % Évaluation perturbée
    [f_xk, c_xk] = fc(xk);
    n_fonc = n_fonc+1;

    % Différences finies avant
    gf(i) = (f_xk-fk)/hi;
    gc(:,i) = (c_xk-ck)/hi;
end
end
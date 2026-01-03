function [gf, gc, n_fonc] = gradient(fc, x, h, nfonc, fx, cx)
% Calcul gradient par differences finies

fk = fx;
ck = cx;
n_fonc = nfonc;
n = length(x);
m = length(ck);
gf = zeros(n,1);
gc = zeros(m,n);

for i=1:n
    xk = x;
    if isscalar(h), hi = h; else, hi = h(i); end
    xk(i) = x(i) + hi;

    [f_xk, c_xk] = fc(xk);
    n_fonc = n_fonc+1;

    gf(i) = (f_xk-fk)/hi;
    gc(:,i) = (c_xk-ck)/hi;
end
end
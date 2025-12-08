function [m, n_fonc] = Merite(fc, x, rho, nfonc)
[fx, cx] = fc(x);
m = fx + rho*norm(cx, 1);
n_fonc = nfonc+1;
end
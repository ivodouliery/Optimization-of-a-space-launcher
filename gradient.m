function [gf, gc, n_fonc] = gradient(fc, x, h, nfonc)
[fk, ck] = fc(x);
n_fonc = nfonc+1;
n = length(x);
m = length(ck);
gf = zeros(n,1);
gc = zeros(m,n);

for i=1:n
    xk = x;
    xk(i) = x(i)+h;
    
    [f_xk, c_xk] = fc(xk);
    n_fonc = nfonc+1; 

    gf(i) = (f_xk-fk)/h;
    gc(:,i) = (c_xk-ck)/h;
end
end
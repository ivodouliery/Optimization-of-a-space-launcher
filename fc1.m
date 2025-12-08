function [f, c] = fc1(x)

%Suivi du nombre d'appel
persistent nfonc;
if isempty(nfonc)
        nfonc = 0;
end
nfonc = nfonc + 1;

c = zeros(2,1);
f = 3*x(1)^2 + 2*x(2)^3;
c(1) = x(1) + x(2) - 3;
c(2) = x(1)^2 + x(2)^2 - 5;
end
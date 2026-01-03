function [f, c] = f_Ariane1(x)

% Données
mu = 1700;
k = [0.1101; 0.1532; 0.2154];
ve = [2647.2; 2922.4; 4344.3];
dV_req = 11527;

Mi = zeros(3,1);
Mf = zeros(3,1);
dV = zeros(3,1);

M_p1 = mu;

for j = 3:-1:1
    m_ej = x(j);
    m_sj = k(j) * m_ej;

    Mf(j) = M_p1 + m_sj;
    Mi(j) = Mf(j) + m_ej;

    dV(j) = ve(j) * log(Mi(j) / Mf(j));

    M_p1 = Mi(j);
end
f = Mi(1);
c = sum(dV) - dV_req;

end

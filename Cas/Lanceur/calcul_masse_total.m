function [Mi, Ms, M0] = calcul_masse_total(Me, Mu, k)
% CALCUL_MASSE_TOTAL Calcule les masses initiales et structurelles du lanceur
% [Mi, Ms, M0] = calcul_masse_total(Me, Mu, k)
%
% Entrées :
%   Me : Vecteur des masses d'ergols [kg] (3x1)
%   Mu : Masse de la charge utile [kg]
%   k  : Vecteur des indices de structure [-] (3x1)
%
% Sorties :
%   Mi : Vecteur des masses totales au début de chaque phase [kg] (3x1)
%   Ms : Vecteur des masses structurelles [kg] (3x1)
%   M0 : Masse totale au décollage (équivalent à Mi(1)) [kg]

nb_etages = length(Me);
Ms = zeros(nb_etages, 1);
Mi = zeros(nb_etages, 1);

M_above = Mu; % Masse portée par le dernier étage (Charge Utile)

for j = nb_etages:-1:1
    % Masse structurelle de l'étage j
    Ms(j) = k(j) * Me(j);

    % Masse totale au début de la phase j
    % Mi(j) = Masse portée + Masse structure + Masse ergols
    Mi(j) = M_above + Ms(j) + Me(j);

    % La masse "au-dessus" de l'étage j-1 est la masse totale de l'étage j
    M_above = Mi(j);
end

M0 = Mi(1);
end

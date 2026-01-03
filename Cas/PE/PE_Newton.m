function [Me, x, lambda, f_opt] = PE_Newton(Vp, Ve, k, Mu)
% PE_Newton Résout le problème d'étagement par la méthode de Newton (3 étages)
% [Me, x, lambda, f_opt] = PE_Newton(Vp, Ve, k, Mu)
%
% Entrées (Optionnelles) :
%   Vp : Vitesse propulsive cible [m/s]
%   Ve : Vecteur des vitesses d'éjection [m/s]
%   k  : Vecteur des indices structurels [-]
%   Mu : Masse de la charge utile [kg]
%
% Sorties :
%   Me     : Vecteur colonne [Me1; Me2; Me3] des masses d'ergols [kg]
%   x      : Vecteur colonne [x1; x2; x3] des ratios de masse
%   lambda : Multiplicateur de Lagrange associé à la contrainte de vitesse
%   f_opt  : Valeur de la fonction objectif (opposé du ratio de charge utile)

% --- 1. Paramètres par défaut ---
if nargin < 1, Vp = 10000; end
if nargin < 2, Ve = [2600, 3000, 4400]; end
if nargin < 3, k  = [0.10, 0.15, 0.20]; end
if nargin < 4, Mu = 1000; end

Omega = k ./ (1 + k);

% --- 2. Algorithme de Newton ---
x3 = 2.5;       % Initialisation
tol = 1e-8;
max_iter = 50;

for i = 1:max_iter
    % Calcul en cascade (x2 et x1 dépendent de x3)
    x2 = (1/Omega(2)) * (1 - (Ve(3)/Ve(2)) * (1 - Omega(3)*x3));
    x1 = (1/Omega(1)) * (1 - (Ve(2)/Ve(1)) * (1 - Omega(2)*x2));

    % Sécurité divergence
    if x1 <= 1e-3 || x2 <= 1e-3
        error('PE_Newton:Divergence', 'Newton a divergé vers des valeurs non physiques (x <= 0).');
    end

    % Évaluation de la contrainte g(x3) = Sum(DeltaV) - Vp
    g = Ve(1)*log(x1) + Ve(2)*log(x2) + Ve(3)*log(x3) - Vp;

    % Convergence ?
    if abs(g) < tol
        break;
    end

    % Dérivée g'(x3)
    term1 = Omega(3) / (Omega(1) * x1);
    term2 = Omega(3) / (Omega(2) * x2);
    term3 = 1.0 / x3;
    gp = Ve(3) * (term1 + term2 + term3);

    % Mise à jour Newton
    x3 = x3 - g / gp;
end

if i == max_iter
    warning('Newton n''a pas convergé après %d itérations.', max_iter);
end

% --- 3. Construction des résultats ---
x = [x1; x2; x3];

% Calcul de f_opt (Opposé du Ratio J)
y3 = (1 + k(3))/x3 - k(3);
y2 = (1 + k(2))/x2 - k(2);
y1 = (1 + k(1))/x1 - k(1);
f_opt = -1 * (y1 * y2 * y3);

% Calcul de Lambda (via étage 1)
C1 = Ve(1) * (1 - Omega(1) * x1);
lambda = f_opt / C1;

% --- 4. Calcul des Masses d'Ergols ---
Me = zeros(3,1);
M_above = Mu;

for j = 3:-1:1
    % Formule dérivée de x_j = (M_above + (1+k)*Me) / (M_above + k*Me)
    Me(j) = M_above * (x(j) - 1) / (1 + k(j) - x(j) * k(j));

    % Mise à jour masse 'au-dessus' pour l'étage suivant (en descendant)
    Ms_j = k(j) * Me(j);
    M_above = M_above + Ms_j + Me(j);
end

M_total = sum(Me) + sum(k'.*Me) + Mu;
J = Mu / M_total;
fprintf('\n--- Résultats PE_Newton ---\n');
fprintf('Masses d''ergols (kg) :\n');
fprintf('Masse totale M0 = %.4f kg\n', M_total);
fprintf('Ratio de Charge Utile calculé : J = %.6f\n', J);

end
% =========================================================================
% RESOLUTION DE L'EQUATION DE CONTRAINTE PAR NEWTON (Section 1.4.4)
% Equation : g(x3) = sum(v_i * ln(x_i)) - Vp = 0
% =========================================================================
clear; clc;

%% 1. PARAMÈTRES
% -------------------------------------------------------------------------
% Vitesse propulsive cible (m/s) - À ajuster selon votre itération globale
Vp = 9341.2;

% Données Lanceur (Ariane 1)
Ve = [2600, 3000, 4400];      % Vitesses d'éjection
k  = [0.10, 0.15, 0.20];      % Indices constructifs
Omega = k ./ (1 + k);         % Paramètre structurel

%% 2. ALGORITHME DE NEWTON
% -------------------------------------------------------------------------
% Initialisation
x3 = 2.5;       % Valeur de départ (doit être > 1)
tol = 1e-8;      % Tolérance de convergence
max_iter = 50;   % Sécurité boucle infinie

fprintf('Début de la méthode de Newton (Cible Vp = %.2f m/s)...\n', Vp);
fprintf('Iter |      x3      |      x2      |      x1      |  Erreur (g) \n');
fprintf('--------------------------------------------------------------\n');

for i = 1:max_iter

    % --- ÉTAPE A : Calcul en cascade (Récurrence) ---
    % On calcule x2 et x1 qui dépendent de x3 actuel
    % x_j = (1/Omega_j) * [ 1 - (Ve_next/Ve_curr)*(1 - Omega_next*x_next) ]

    % x2 en fonction de x3
    x2 = (1/Omega(2)) * (1 - (Ve(3)/Ve(2)) * (1 - Omega(3)*x3));

    % x1 en fonction de x2
    x1 = (1/Omega(1)) * (1 - (Ve(2)/Ve(1)) * (1 - Omega(2)*x2));

    % Sécurité : Logarithme impossible si x <= 0
    if x1 <= 1e-3 || x2 <= 1e-3
        error('Newton a divergé vers des valeurs non physiques (x <= 0). Changez x3 initial.');
    end

    % --- ÉTAPE B : L'Équation (La Contrainte g(x3)) ---
    % C'est l'équation que l'on veut annuler : Somme(DeltaV) - Vp = 0
    g = Ve(1)*log(x1) + Ve(2)*log(x2) + Ve(3)*log(x3) - Vp;

    % --- ÉTAPE C : La Dérivée g'(x3) ---
    % Formule analytique calculée précédemment
    term1 = Omega(3) / (Omega(1) * x1);
    term2 = Omega(3) / (Omega(2) * x2);
    term3 = 1.0 / x3;

    gp = Ve(3) * (term1 + term2 + term3);

    % Affichage suivi
    fprintf('%4d | %10.6f | %10.6f | %10.6f | %10.2e\n', i, x3, x2, x1, g);

    % --- ÉTAPE D : Mise à jour et Test ---
    if abs(g) < tol
        fprintf('--------------------------------------------------------------\n');
        fprintf('>>> CONVERGENCE ATTEINTE en %d itérations.\n', i);
        break;
    end

    % Formule de Newton : x_new = x_old - g(x) / g'(x)
    x3 = x3 - g / gp;
end

%% 3. RÉSULTATS FINAUX
% -------------------------------------------------------------------------
fprintf('\nSolution optimale trouvée :\n');
fprintf('x1 = %.6f\n', x1);
fprintf('x2 = %.6f\n', x2);
fprintf('x3 = %.6f\n', x3);

% Vérification finale
V_reelle = Ve(1)*log(x1) + Ve(2)*log(x2) + Ve(3)*log(x3);
fprintf('Vitesse atteinte = %.4f m/s (Cible : %.4f)\n', V_reelle, Vp);

%% 4. VÉRIFICATION DES CONDITIONS KKT ET CALCUL DE LAMBDA
% -------------------------------------------------------------------------
fprintf('\n--- VERIFICATION KKT & MULTIPLICATEUR ---\n');

% 1. Calcul de la fonction objectif f au point optimum
% f = - y1 * y2 * y3 (Opposé du ratio de performance global)
y3 = (1 + k(3))/x3 - k(3);
y2 = (1 + k(2))/x2 - k(2);
y1 = (1 + k(1))/x1 - k(1);

f_opt = -1 * (y1 * y2 * y3);
fprintf('Valeur du critère f(x*) : %.6f\n', f_opt);

% 2. Vérification de la constante structurelle pour chaque étage
% C_j = Ve_j * (1 - Omega_j * x_j)
% D'après la théorie, C_1 doit être égal à C_2 et C_3.

C = zeros(1, 3); % Vecteur pour stocker les constantes
Lambda = zeros(1, 3); % Vecteur pour stocker les lambda calculés

for j = 1:3
    % Calcul de la constante structurelle
    C(j) = Ve(j) * (1 - Omega(j) * eval(sprintf('x%d', j)));

    % Calcul du multiplicateur lambda associé
    % Relation : Lambda * C_j = f_opt  =>  Lambda = f_opt / C_j
    Lambda(j) = f_opt / C(j);
end

% 3. Affichage des résultats
fprintf('\nVérification de la constante v_ej * (1 - Omega_j * x_j) :\n');
fprintf('Etage 1 : %.4f\n', C(1));
fprintf('Etage 2 : %.4f\n', C(2));
fprintf('Etage 3 : %.4f\n', C(3));

fprintf('\nCalcul du multiplicateur de Lagrange (Lambda) :\n');
fprintf('Lambda (via Etage 1) : %.4e\n', Lambda(1));
fprintf('Lambda (via Etage 2) : %.4e\n', Lambda(2));
fprintf('Lambda (via Etage 3) : %.4e\n', Lambda(3));

% Conclusion automatique
ecart_relatif = std(Lambda) / mean(abs(Lambda));
if ecart_relatif < 1e-4
    fprintf('\n>>> SUCCÈS : Les conditions KKT sont vérifiées.\n');
    fprintf('    Le multiplicateur est unique (à la précision près).\n');
    fprintf('    Valeur retenue : Lambda = %.4e\n', mean(Lambda));
else
    fprintf('\n>>> ATTENTION : Les conditions KKT ne semblent pas vérifiées.\n');
end
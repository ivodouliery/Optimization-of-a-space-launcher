clear;
close all;

% === CONSTANTES PHYSIQUES ===
mu_terre = 3.986e14;  % [m^3/s^2] Constante gravitationnelle
cx      = 0.1;       % [-] Coefficient de trainee
rho0    = 1.225;     % [kg/m^3] Densite atmospherique au sol
H       = 7000;      % [m] Hauteur d'echelle
R_terre  = 6378137;   % [m] Rayon terrestre
g0      = 9.81;      % [m/s^2] (Si besoin)
Altitude_cible = 200000; % [m] Altitude cible
Rc      = R_terre + Altitude_cible;   % [m] Rayon orbite cible

% === PARAMETRES LANCEUR ===
% Charge Utile
Mu = 1000; % [kg]

% Parametres propulsifs des etages (1 -> 2 -> 3)
Alpha = [15; 10; 10];      % [m/s^2] Acceleration initiale demandee
Vc    = sqrt(mu_terre/Rc); % [m/s] Vitesse orbitale cible
Vp    = 9251;             % [m/s] Vitesse propulsive cible
Ve    = [2600; 3000; 4400]; % [m/s] Vitesse d'ejection
k     = [0.1; 0.15; 0.2];  % [-] Indices de structure (Hypothese ou valeurs Ariane 1)

% === ETAGEAMENT OPTIMAL ===
% Calcul des masses d'ergols (Me) et ratios via Newton
% Me est directement retourne par la fonction modifiee
Me = PE_Newton(Vp, Ve, k, Mu);


% === CALCUL DES PERFORMANCES ===
% Initialisation
T  = zeros(3,1); % Poussee [N]
q  = zeros(3,1); % Debit massique [kg/s]
tc = zeros(3,1); % Temps de combustion [s]

M_above = Mu; % Masse portee par le dernier etage

% Calcul des masses structurelles et initiales
[Mi, Ms, M0] = calcul_masse_total(Me, Mu, k);

for j = 1:3
    % Calcul Poussee (T = alpha * Mi)
    T(j) = Alpha(j) * Mi(j);

    % Calcul Debit (q = T / Ve)
    q(j) = T(j) / Ve(j);

    % Duree de combustion (tc = Me / q)
    tc(j) = Me(j) / q(j);
end

x0 = [0.02; 0.05; 0.05; 0.05];
lambda0 = [0; 0];
ub = [0.0;  0.3;  0.5;  0.5];
lb = [0.0 ; 0.3; -0.5; -0.5];

options.hessian_mod = 1;
options.max_iter = 100;
options.max_eval = 500;
options.tol_x = 1e-8;
options.tol_f = 1e-8;

param = struct('mu_terre', 3.986e14, 'cx', 0.1, 'rho0', 1.225, 'H', 7000, 'R_terre', 6378137, 'g0', 9.81, 'Mi', Mi, 'Ms', Ms, 'M0', M0, 'tc', tc, 'T', T, 'q', q, 'Rc', Rc);

% Facteurs d'echelle pour l'optimisation
% On veut ramener f et c a l'ordre de grandeur de 1
scale_f = 10000;      % Vitesse ~ 7500 m/s -> 7.5
scale_c1 = 10000;    % Distance ~ DeltaR -> 10^4 m
scale_c2 = 1e8;      % R.V ~ 6e6 * 100 -> 6e8

fc = @(x) f_wrapper(x, param, scale_f, scale_c1, scale_c2);

% BALAYAGE (Theta0, Theta1) ---
x0 = balayage(fc, x0);

% OPTIMISATION PARTIELLE (Theta0, Theta1 fixés) ---
lb_fixed = lb;
ub_fixed = ub;
% On bloque les deux premiers parametres sur la valeur trouvee
lb_fixed(1:2) = x0(1:2);
ub_fixed(1:2) = x0(1:2);

% On lance l'SQP sur les autres variables
[x_found_1, lambda_1] = SQP(fc, x0, lambda0, lb_fixed, ub_fixed, options);

% RAFFINEMENT GLOBAL ---
x0_phase2 = x_found_1;
lambda0_phase2 = lambda_1;

% On remet les bornes originales (avec theta0 >= 0 strict)
% lb doit etre coherent avec ce qu'on veut
lb_final = [0.0; -0.6; -0.5; -0.5];
ub_final = [pi/2; 0.6; 0.5; 0.5];

[theta, lambda] = SQP(fc, x0_phase2, lambda0_phase2, lb_final, ub_final, options);

% Wrapper pour mise a l'echelle
function [f_sc, c_sc] = f_wrapper(x, p, sf, sc1, sc2)
[f, c] = f_trajectoire(x, p);
f_sc = f / sf;
c_sc = [c(1)/sc1; c(2)/sc2];
end

R = [param.R_terre; 0];
V = 100*[cos(theta(1)); sin(theta(1))];
M = M0;

opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
tcurr = 0;

% Stockage pour tout le vol
t_cell = cell(3, 1);
y_cell = cell(3, 1);

for j = 1:3
    % Integration de l'etage j
    tspan = [tcurr, tcurr + tc(j)];

    % Attention : theta(1) est l'angle initial, theta(2) est l'incidence etage 1, etc.
    [t, y] = ode45(@(t,y) equation_mouvement(t,y,theta(j+1),T(j),q(j),param), tspan, [R; V; M], opts);

    % Stockage des resultats dans les cellules
    t_cell{j} = t;
    y_cell{j} = y;

    % Mise a jour etat initial pour etage suivant
    % y(end,:) est une ligne, on veut des colonnes pour R et V
    R = y(end,1:2)';
    V = y(end,3:4)';

    % Masse fin de combustion - Masse structurelle = Masse debut etage suivant
    M = y(end,5) - Ms(j);

    tcurr = t(end);
end

% Assemblage final des resultats
t_all = vertcat(t_cell{:});
y_all = vertcat(y_cell{:});

% === AFFICHAGE DES RESULTATS ===

% Figure 1 : Trajectoire SCHEMATIQUE (Normalisée)
% Le but est de bien voir la trajectoire par rapport a la Terre
% sans que la Terre ecrase tout (car 200km << 6400km).
figure(1); clf; hold on; axis equal; grid on;

% Parametres graphiques
R_terre_plot = 1.0;
H_orbite_plot = 0.3; % On represente l'orbite a 30% du rayon terrestre pour la visibilite
Scale_Alt = H_orbite_plot / Altitude_cible; % Facteur de dilatation de l'altitude

% Transformation des donnees physiques -> Données de plot
% 1. On passe en polaire (r, theta)
r_phys = sqrt(y_all(:,1).^2 + y_all(:,2).^2);
theta_phys = abs(atan2(y_all(:,2), y_all(:,1))); % Force l'angle positif pour affichage vers le "Haut"

% 2. On applique le scaling sur l'altitude (r - Rt)
alt_phys = r_phys - param.R_terre;
r_plot = R_terre_plot + alt_phys * Scale_Alt;

% 3. On repasse en cartesien pour le plot (X_plot, Y_plot)
% On garde l'orientation: X vertical, Y horizontal (mais plot inverse X/Y apres)
X_plot = r_plot .* cos(theta_phys);
Y_plot = r_plot .* sin(theta_phys);

% --- AFFICHAGE ---

% 1. Terre (Disque unitaire)
theta_circle = linspace(0, 2*pi, 200);
fill(R_terre_plot*cos(theta_circle), R_terre_plot*sin(theta_circle), [0.85 0.85 0.85], 'EdgeColor', 'k', 'LineWidth', 1);
text(0, 0, 'Terre', 'HorizontalAlignment', 'center', 'FontSize', 12);

% 2. Orbite Cible (Cercle ray 1.3)
R_c_plot = R_terre_plot + H_orbite_plot;
plot(R_c_plot*cos(theta_circle), R_c_plot*sin(theta_circle), 'k--', 'LineWidth', 1.5);
text(0, R_c_plot + 0.05, 'Orbite Cible', 'HorizontalAlignment', 'center');

% 3. Trajectoire
plot(X_plot, Y_plot, 'b-', 'LineWidth', 2); % X_plot en abscisse, Y_plot en ordonnee (Depart a Droite)

% 4. Point final
plot(X_plot(end), Y_plot(end), 'b.', 'MarkerSize', 15);

title('Trajectoire Schématique (Altitude Exagérée)');
xlabel('Altitude X Normalisée'); ylabel('Portée Y Normalisée');
legend('Terre', 'Orbite Cible', 'Trajectoire', 'Location', 'bestoutside');
axis([-1.5 1.5 -1.5 1.5]); % Vue d'ensemble fixe

% Figure 2 : Altitude en fonction du temps
figure(2); clf; hold on; grid on;
altitude_h = sqrt(y_all(:,1).^2 + y_all(:,2).^2) - param.R_terre;
plot(t_all, altitude_h, 'b', 'LineWidth', 1.5);
yline(Altitude_cible, 'r--', 'Altitude Cible');
title('Altitude vs Temps');
xlabel('Temps [s]'); ylabel('Altitude [m]');

% Figure 3 : Vitesse en fonction du temps
figure(3); clf; hold on; grid on;
vitesse_abs = sqrt(y_all(:,3).^2 + y_all(:,4).^2);
plot(t_all, vitesse_abs, 'm', 'LineWidth', 1.5);
yline(Vc, 'r--', 'Vitesse Cible');
title('Vitesse vs Temps');
xlabel('Temps [s]'); ylabel('Vitesse [m/s]');

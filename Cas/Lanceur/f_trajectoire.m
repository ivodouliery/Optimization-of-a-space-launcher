function [f, c] = f_trajectoire(x, param)
% f_trajectoire Calcule le critère et les contraintes pour l'optimisation de la trajectoire
% [f, c] = f_trajectoire(x, param)
%
% Entrées :
%   x     : Vecteur des variables d'optimisation (4x1)
%           x(1) = theta_0 (Angle vitesse initiale) [rad]
%           x(2) = theta_1 (Angle poussée Etage 1) [rad]
%           x(3) = theta_2 (Angle poussée Etage 2) [rad]
%           x(4) = theta_3 (Angle poussée Etage 3) [rad]
%   param : Structure contenant les paramètres du lanceur et physiques
%           .mu_terre, .R_terre, .cx, .rho0, .H
%           .T (3x1), .q (3x1), .tc (3x1)
%           .Mi (3x1), .Ms (3x1)
%           .Rc (Rayon orbite cible)
%
% Sorties :
%   f : Coût à minimiser (-V_finale à l'injection)
%   c : Vecteur des contraintes d'égalité [c_orbite; c_rdv]
%       c(1) = R_final - Rc
%       c(2) = Dot(R_final, V_final) (Condition d'injection horizontale)

% --- 1. Initialisation ---
theta0 = x(1);

% Position initiale : Surface terrestre sur l'axe X
r0 = [param.R_terre; 0];

% Vitesse initiale : V0 = 100 m/s orientée par theta0
% (Hypothèse "lift-off" pour éviter la singularité à V=0)
v0 = 100 * [cos(theta0); sin(theta0)];

% Masse initiale (Au décollage)
m0 = param.Mi(1);

% Etat initial global [rx; ry; vx; vy; m]
y0 = [r0; v0; m0];

t_curr = 0;
y_curr = y0;

% --- 2. Simulation des 3 étages ---
for j = 1:3
    % Paramètres de l'étage j
    T_j = param.T(j);
    q_j = param.q(j);
    tc_j = param.tc(j);
    theta_j = x(j+1); % theta_1, theta_2, theta_3

    % Tolerances adaptees aux composants [rx ry vx vy m]
    opts = odeset('RelTol', 1e-6, 'AbsTol', [1; 1; 1e-3; 1e-3; 0.1]);

    [t_out, Y_stage] = ode45(@(t, y) equation_mouvement(t, y, theta_j, T_j, q_j, param), ...
        [t_curr, t_curr + tc_j], y_curr, opts);

    % Etat à la fin de la combustion
    y_curr = Y_stage(end, :)';
    t_curr = t_out(end);

    if j < 3
        y_curr(5) = param.Mi(j+1);
    end
end

% --- 3. Calcul Coût et Contraintes ---
r_f = y_curr(1:2);
v_f = y_curr(3:4);

R_final = norm(r_f);
V_final = norm(v_f);

% DEBUG:
% fprintf('Resultat Final : V_final = %.4f m/s\n', V_final);

% Critère : Maximiser V_finale <=> Minimiser -V_finale
f = -V_final;

% Contraintes
% 1. Altitude orbitale
c1 = R_final - param.Rc;

% 2. Injection horizontale (R orthogonal à V)
c2 = dot(r_f, v_f);

c = [c1; c2];

end

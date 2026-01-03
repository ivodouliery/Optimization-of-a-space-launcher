function dy = equation_mouvement(~, y, theta_incidence, T, q, param)
% y = [rx; ry; vx; vy; m]
rx = y(1);
ry = y(2);
vx = y(3);
vy = y(4);
m  = y(5);

r_norm = sqrt(rx^2 + ry^2);
v_norm = sqrt(vx^2 + vy^2);

% Gravité
gx = -param.mu_terre * rx / r_norm^3;
gy = -param.mu_terre * ry / r_norm^3;

% Traînée (Modèle du PDF)
alt = r_norm - param.R_terre;
if alt < 0, rho = param.rho0; else, rho = param.rho0 * exp(-alt / param.H); end

Dx = -param.cx * rho * v_norm * vx;
Dy = -param.cx * rho * v_norm * vy;

% Poussée
% Calcul de la pente de la vitesse (gamma)
gamma = atan2(vy, vx);

% L'angle absolu de poussée est la somme de la pente vitesse et de l'incidence
angle_poussee = gamma + theta_incidence;

Tx = T * cos(angle_poussee);
Ty = T * sin(angle_poussee);

% Somme des forces
ax = (Tx + Dx) / m + gx;
ay = (Ty + Dy) / m + gy;

dy = [vx; vy; ax; ay; -q];
end
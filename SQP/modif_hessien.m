function [Q_prime] = modif_hessien(Q)
% Régularisation du Hessien (si non défini positif)

lambda_min = min(eig(Q));
if lambda_min <= 0
    % Ajout d'un terme diagonal pour rendre Q défini positif
    tho = -lambda_min + 1e-6;
else
    tho = 0;
end

Q_prime = Q + tho * eye(size(Q,1));
end
function [Q_prime] = modif_hessien(Q)
% Regularisation Hessien (si valeurs propres negatives)

lambda_min = min(eig(Q));
if lambda_min <= 0
    tho = -lambda_min + 1e-6;
else
    tho = 0;
end

Q_prime = Q + tho * eye(size(Q,1));
end
function [lambda_QP, d_QP] = Sol_pb_quadra(Q, gf, gc, ck)
% Résolution du sous-problème quadratique (QP)

g = gf;
A = gc;
b = -ck;

% Résolution du système linéaire (pinv pour robustesse face aux singularités)
M = A * (Q \ A');
lambda_QP = - pinv(M) * (A * (Q \ g) + b);

% Direction de descente
d_QP = -Q\(A'*lambda_QP + g);

end
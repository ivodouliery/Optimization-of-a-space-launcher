function [lambda_QP, d_QP] = Sol_pb_quadra(Q, gf, gc, ck)
% Resolution du sous-probleme quadratique (QP)

g = gf;
A = gc;
b = -ck;

% Resolution (pinv pour robustesse)
M = A * (Q \ A');
lambda_QP = - pinv(M) * (A * (Q \ g) + b);
d_QP = -Q\(A'*lambda_QP + g);

end
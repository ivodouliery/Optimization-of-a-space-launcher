function [lambda_QP, d_QP] = Sol_pb_quadra(Q, gf, gc, ck)
g = gf;
A = gc;
b = -ck;

lambda_QP = - (A * (Q \ A')) \ (A * (Q \ g) + b);
d_QP = -Q\(A'*lambda_QP + g);

end
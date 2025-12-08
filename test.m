clear;
close all;


x0 = [-1;2;1;-2;-2];
lambda0 = zeros(3,1);
h = 1e-10;

[x_star, lambda_star] = SQP(@f_MHW4D, x0, lambda0, h, 1);
disp(x_star);




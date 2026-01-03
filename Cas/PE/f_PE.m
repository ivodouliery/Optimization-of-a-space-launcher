function [f, c] = f_PE(x, k, Ve, Vp)

term1 = (1 + k(1)) / x(1) - k(1);
term2 = (1 + k(2)) / x(2) - k(2);
term3 = (1 + k(3)) / x(3) - k(3);

f = -(term1 * term2 * term3);

c = (Ve(1)*log(x(1)) + Ve(2)*log(x(2)) + Ve(3)*log(x(3)) - Vp) / 30000;

end
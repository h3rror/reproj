function r_star_ = r_star(r) % only correct for linear + quadratic term

r_star_ = r + (r+1).*r/2;
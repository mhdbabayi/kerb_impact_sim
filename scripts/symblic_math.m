syms rho alpha gamma D A B beta phi real
clc
eps = exp(-beta*phi)*(A*sin(beta*phi) + B*cos(beta*phi));
normal_force = (int(eps*cos(phi)));
eps_d1 = diff(eps, phi);
eps_d2 = diff(eps_d1 , phi);
eps_exp  = eval(subs(eps_d2 , phi , 0));
alpha = asin((D/rho)*sin(gamma))- gamma;
Rp = rho*sin(alpha)/sin(gamma);
d1 = diff(Rp , gamma);
d2 = diff(d1 , gamma);
d2_simp = subs(d2 , (gamma - asin(D*sin(gamma)/rho)), 'alpha');
contact_patch_integral = int(A * phi^2*cos(phi));
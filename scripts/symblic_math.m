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
%%
clc
d5 = diff(-Rp, gamma , 5);
d5_simp = simplify(d5);
d5_simp_sub = subs(d5_simp , sin(gamma) , "S");
d5_simp_sub = subs(d5_simp_sub, cos(gamma) ,1)
d5_f = matlabFunction(subs(d5_simp_sub, {"rho", "D"} , {0.5, 0.7}))
%%
clear 
clc
close all
syms rho_s D_s R_s kt real positive
syms phi lambda C0 C1 C2 real
tyre_radius = 0.778/2;
penetration = 0.1;
rho = 0.07;
R0 = tyre_radius + rho - penetration;
limit_theta = acos((R0^2 + tyre_radius^2 - rho^2)/(2*R0*tyre_radius));
full_theta = linspace(0 , 2*pi, 360);
theta = linspace(-limit_theta , limit_theta, 100);
terrain_x = R0 + rho*cos(full_theta);
terrain_y = rho*sin(full_theta);
polar_circle = D_s*cos(phi) - sqrt(D_s^2*cos(phi)^2 +(rho_s)^2 - D_s^2);
circle_taylor = taylor(polar_circle, phi, "order", 6);
RHS = -diff(circle_taylor, phi , 5) - 2*diff(circle_taylor, phi, 3) - (lambda + kt)*diff(circle_taylor, phi);
% The coefficients from RHS:
D = -(D_s/24 - D_s^2/(6*rho) + D_s/(8*rho_s^3));
E1 = 4*D*(kt + lambda);
E2 = 4*D*(kt + lambda)+((kt+lambda)*(kt+lambda)*(D_s/2 - D_s^2/(2*rho))) + 48*D;
A = E1/(4*lambda);
B = E2/(2*lambda) - 12*E1/(2*lambda^2);
alpha = sqrt((lambda-1)/2);
beta = sqrt((lambda + 1)/2);
t = C1*cosh(alpha*phi)*cos(beta*phi) + ...
    C2*sinh(alpha*phi)*sin(beta*phi) + ...
    A * phi^4 + B*phi^2 + C0;
R_s = tyre_radius - (circle_taylor + t);
belt_deflection = circle_taylor+ t;
bc_1 = subs(diff(belt_deflection , phi , 2), phi , 0);
final_radius_function = matlabFunction(subs(R_s,[C0, C1, C2, D_s, kt, lambda, rho_s],...
                                  {0.01, 0.1, 0.1, tyre_radius + 10 - 0.1, 100000, 5000, 10}));
close all
smaple_phi = linspace(-pi/6, pi/6, 61);
plot(rad2deg(smaple_phi), final_radius_function(smaple_phi))
%%
plot(terrain_x , terrain_y);
hold on
plot(tyre_radius*cos(full_theta), tyre_radius*sin(full_theta) , 'bla')
circle_taylor_f = matlabFunction(subs(circle_taylor, {D_s , rho_s}, {R0, rho}));
R_taylor = circle_taylor_f(theta);
hold on
plot(R_taylor.*cos(theta) , R_taylor.*sin(theta), '*-')
daspect([1 1 1])

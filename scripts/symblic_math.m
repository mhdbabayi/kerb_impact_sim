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
d5_simp_sub = subs(d5_simp_sub, cos(gamma) ,1);
d5_f = matlabFunction(subs(d5_simp_sub, {"rho", "D"} , {0.5, 0.7}));
%%
clear 
clc
close all
syms rho Dc R kt real positive
syms phi lambda lambdaT C0 C1 C2 real
constants = struct ('R' , 0.778/2, 'penetration', 0.1,...
    'rho', 1, 'kt', 1.4567e8,'EI', 13.9,'kr',1.5e6,...
    'C0', 0.01, 'C1', 0.000001, 'C2', 0.0001);

constants.D = constants.R + constants.rho - constants.penetration;
constants.lambda = sqrt(constants.kr*constants.R^4/constants.EI + 1);
constants.lambdaT = sqrt(constants.kt*constants.R^4/constants.EI);
    limit_theta = acos((constants.D^2 + constants.R^2 - constants.rho^2)/(2*constants.D*constants.R));
constants.phi = linspace(-limit_theta , limit_theta , 61);
my_eval = @(exp)vpa(subs(exp , [R, Dc,lambda ,lambdaT , rho, C0, C1 , C2],...
    {constants.R, constants.D, constants.lambda , constants.lambdaT, constants.rho,...
    constants.C0 , constants.C1, constants.C2}));

polar_circle = Dc*cos(phi) - sqrt(Dc^2*cos(phi)^2 +(rho)^2 - Dc^2);
circle_taylor = taylor(polar_circle, phi, "order", 6);
RHS = diff(circle_taylor, phi , 5) + 2*diff(circle_taylor, phi, 3) + (lambda^2)*diff(circle_taylor, phi);
% The coefficients from RHS:
D = (Dc/24 - Dc^2/(6*rho) + Dc^4/(8*rho^3));
E1 = 4*D*(lambda^2);
E2 = (lambda^2) * ((Dc^2/rho - Dc)) + 48*D;
A = D * (lambda^2)/(lambda^2 + lambdaT^2);
B = E2/(lambda^2 + lambdaT^2) - 24*D/(lambda^2 + lambdaT)^2;
alpha = sqrt((lambda-1)/2);
beta = sqrt((lambda + 1)/2);
general_solution = C1*cosh(alpha*phi)*cos(beta*phi) + ...
                    C2*sinh(alpha*phi)*sin(beta*phi);
particular_solution = A * phi^4 + B*phi^2 + C0;
diff_eq_operator = @(sol)diff(sol , phi , 5) + 2*diff(sol , phi , 3) + lambda*(diff(sol , phi));
t = general_solution + particular_solution;
final_radius = circle_taylor -t;
belt_deflection = circle_taylor -t;
bc_1 = subs(diff(belt_deflection , phi , 2), phi , 0);
T_f = matlabFunction(eval(subs(t ,[C0 , C1 , C2 , Dc , lambda, lambdaT, rho],...
    {constants.C0 , constants.C1 , constants.C2 , constants.D, constants.lambda, constants.lambdaT , constants.rho})));
final_radius_function = matlabFunction(subs(final_radius,[C0, C1, C2, Dc, lambdaT, lambda, rho],...
                                 {constants.C0,constants.C1, constants.C2, constants.D, constants.lambdaT, constants.lambda, constants.rho}));
close all
set(0 , 'DefaultLineLineWidth' , 2)
theta = linspace(-pi , pi , 361);
terrain_x = constants.rho * cos(theta);
terrain_y = constants.rho*sin(theta) - constants.D;
plot(terrain_x, terrain_y);
tyre_x = constants.R*cos(theta);
tyre_y = constants.R*sin(theta);
hold on
plot(tyre_x, tyre_y);
daspect([1 1 1])
fitted_circle_R = eval(subs(circle_taylor , [Dc, rho , phi], {constants.D , constants.rho , constants.phi}));
fitted_circle_y = -fitted_circle_R .* cos(constants.phi);
fitted_circle_x = fitted_circle_R .* sin(constants.phi);
plot(fitted_circle_x , fitted_circle_y, LineWidth=3, Color='magenta')
deflected_R = final_radius_function(constants.phi);
deflected_x = deflected_R .* sin(constants.phi);
deflected_y = -deflected_R .* cos(constants.phi);               
plot(deflected_x , deflected_y)
%%
plot(terrain_x , terrain_y);
hold on
plot(tyre_radius*cos(full_theta), tyre_radius*sin(full_theta) , 'bla')
circle_taylor_f = matlabFunction(subs(circle_taylor, {Dc , rho}, {R0, rho}));
R_taylor = circle_taylor_f(theta);
hold on
plot(R_taylor.*cos(theta) , R_taylor.*sin(theta), '*-')
daspect([1 1 1])

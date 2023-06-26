close all
clc
clear
% R : tyre_radius, rho = terrain radius
%% constants
R= 0.788/2;
R_rim = 0.508/2;
rho = R*0.5;
penetration = 0.03;
D = R + rho - penetration;
theta0 = deg2rad(-45);
base_theta = deg2rad((0:360));
Ot = [D*cos(theta0) , D*sin(theta0)];
max_theta_p = acosd(-(rho^2 - (R^2 + D^2))/(2*R*D));
beta = 1.3;
%% symbolics
clc
syms theta_s A_s B_s phi theta_p
theta_f = deg2rad((0.1:0.1:max_theta_p));
alpha = asin(D*sin(theta_p)/rho) - theta_p;
epsilon0_s = R - sin(alpha)*rho/sin(theta_p);
epsilon0_f = matlabFunction(epsilon0_s);
d1_epsilon0_s = diff(epsilon0_s , theta_p);
d1_epsilon0_f = matlabFunction(d1_epsilon0_s);
d2_epsilon_s = diff(d1_epsilon0_s , theta_p);
d2_epsilon0_f = matlabFunction(d2_epsilon_s);
beam_s = exp(-beta*phi)*(A_s*sin(beta*phi) + B_s*cos(beta*phi));
d1_beam_s = diff(beam_s , phi);
d2_beam_s = diff(d1_beam_s, phi);
d2_beam_f = matlabFunction(subs(d2_beam_s, phi , 0));
epsilon0 = R - epsilon0_f(theta_f);
A_f = d1_epsilon0_f(theta_f)/beta + epsilon0;
beam_d2 = d2_beam_f(A_f);
epsilon0_d2 = d2_epsilon0_f(theta_f);
plot(rad2deg(theta_f) , beam_d2);
hold on
plot(rad2deg(theta_f) , epsilon0_d2)
separation_angle = theta_f(find(beam_d2 > epsilon0_d2 , 1 , 'first'));
if isempty(separation_angle) || separation_angle > rad2deg(max_theta_p)
    separation_angle = deg2rad(max_theta_p);
end
R_separation = R - epsilon0_f(separation_angle);
rad2deg(separation_angle)
beam_theta = deg2rad((0:1:120));
beam_f = matlabFunction(beam_s);
beam_profile_f = @(theta) beam_f(d1_epsilon0_f(separation_angle)/beta + epsilon0_f(separation_angle),...
    epsilon0_f(separation_angle), theta);
beam_R = R -  beam_profile_f(beam_theta);
legend beam terrain
%% graphics objects denoted with _l
figure
hold on
set(0 , "DefaultLineLineWidth", 3)
set(0 , "DefaultAxesFontSize", 20)
set(gca ,"fontsize", 30)
tyre_l = line(R*cos(base_theta) , R*sin(base_theta) , "color" , "b", "linestyle", "--");
terrain_l = line(Ot(1)+ rho*cos(base_theta), Ot(2) + rho*sin(base_theta), "color", "black");
rim_l = line(R_rim*cos(base_theta) , R_rim*sin(base_theta), "linestyle", "-.");
contact_centre_l = line((R- penetration)*cos(theta0) , (R - penetration)*sin(theta0),...
    "color", "black","marker", "*","markersize", 20);
draw_coordinate_frame_polar(0 , 0 , [1 , 0] , 0.1*R, "TCF")
draw_coordinate_frame_polar(Ot(1) , Ot(2) , [sin(theta0), cos(theta0)], 0.1*R , "CCF")
separation_point = line(R_separation*cos(theta0 + separation_angle) , R_separation*sin(theta0+ separation_angle) , 'Marker',"*");
beam_l = line(beam_R.*cos(theta0+ separation_angle + beam_theta) , beam_R.*sin(theta0+separation_angle+beam_theta),...
    "linestyle", ":", "color", "red", "marker", ".", "markersize", 3);
original_springs_l = [];
% for t = 270:10:360
%     original_springs_l(end+1) = line([tyre_l.XData(t), rim_l.XData(t)], [tyre_l.YData(t), rim_l.YData(t)]);
% end
daspect([1, 1, 1])
%%
function draw_coordinate_frame_cartisian(centre_x, centre_y, x_direction, scale , name)
x_direction = x_direction/(norm(x_direction));
plot(centre_x , centre_y , "bla.", MarkerSize=20);
quiver(centre_x , centre_y , x_direction(1) , x_direction(2),scale, "linewidth" , 3 , "color", "black","maxheadsize", 5);
quiver(centre_x , centre_y , -x_direction(2) , x_direction(1) , scale, "linewidth" , 3, "color", "black", "maxheadsize", 5);
text(centre_x - scale*x_direction(1) , centre_y - scale*x_direction(2), name);
text(centre_x + 1.1*scale*x_direction(1) , centre_y + 1.1*scale*x_direction(2) , "X")
text(centre_x - 1.1*scale*x_direction(2) , centre_y + 1.1*scale*x_direction(1) , "Y")
end
function draw_coordinate_frame_polar(centre_x , centre_y , x_direction, scale, name)
x_direction = x_direction/(norm(x_direction));
plot(centre_x , centre_y , "bla.", MarkerSize=20);
quiver(centre_x , centre_y , x_direction(1) , x_direction(2),scale, "linewidth" , 3 , "color", "black","maxheadsize", 5);
quiver(centre_x , centre_y , -x_direction(2) , x_direction(1) , scale, "linewidth" , 3, "color", "black", "maxheadsize", 5);
text(centre_x - 2*scale*x_direction(1) , centre_y, name, FontSize=18);
text(centre_x + 1.1*scale*x_direction(1) , centre_y + 1.1*scale*x_direction(2) , "\theta_0", FontSize=18)
end

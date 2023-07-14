close all
clc
if (~exist('carData', 'var'))
    load data/kerbImpact_experiment.mat;
end
load data/step_out.mat
%% 
%%% terminology: 
%%% test: a fixed configuration, in this case of speed and pressure.  
%%% run: a single experiment, a test can include multiple runs
clc 
free_radius = (0.788/2); %based on 255/55-R20 tyre
step_height = 0.08; % based on image, not accurate

data_objects  = {};
for i = 1:length(carData(: , 1))
    data_objects{i} = KerbImpactData(carData , i , channel_idx);
end
%%
clc
close all
figure()
fz_linespec = "-b";
spring_height_line_spec = "-r";
for i = 1:length(data_objects)
    yyaxis left
    fz_line_obj = data_objects{i}.quad_plot("fl", "fz", fz_linespec);
    yyaxis right
    spring_height_line_obj = data_objects{i}.quad_plot("fl", "spring_height", spring_height_line_spec);
end
sim_suspension_displacement_cm = (sprung_mass.position.y - tyre.position.y)*100;
sim_suspension_displacement_cm = sim_suspension_displacement_cm - sim_suspension_displacement_cm(1) - 10;
for i = 1:4
    subplot(4 , 1 ,i)
    yyaxis left
    total_force_plot = plot((tyre.position.x - 5), ...
        tyre.force.y+ rigid_ring.force.y, "bla", LineWidth=3, Marker="none");
    yyaxis right
    suspension_force_plot = plot(tyre.position.x-5, ...
        sim_suspension_displacement_cm, 'g--',LineWidth=2 );
end
linkAx()
xlim ([-2*free_radius , 5*free_radius])
legend ([total_force_plot, fz_line_obj, suspension_force_plot, spring_height_line_obj ],...
    'sim_fz', 'experiment_fz', 'sim_suspension', 'experiment_suspension')



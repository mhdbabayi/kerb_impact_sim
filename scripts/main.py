import flexring as flx
import physics_engine as phsx
import time
import os, sys
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from euclid3 import Vector2
os.system("clear")
from scipy import io
# constants and initial values
DEBUG = sys.gettrace()
repo_root_path = Path(__file__).resolve().parent.parent
matlab_file_path = Path.joinpath(repo_root_path, "data", "beta_5.mat") 
output_file_path = Path.joinpath(repo_root_path, "data", "step_out.mat")
forward_speed_kph = 12.
forward_speed = forward_speed_kph/3.6
tyre_radius = 0.788/2
if DEBUG:
    draw_frequency = 1
    initial_x = 4.6
else:
    draw_frequency = 10
    initial_x = 4
sprung_mass = 700
unsprung_mass = 50
initial_y = tyre_radius - (sprung_mass + unsprung_mass)*10/(flx.ContinousTyre.lump_stiffness)
# defining sim objects, all moving objects inherit from rigid body
step_road: flx.Road = flx.Road.make_simple_road(
                step_width=0.02,
                step_height=0.08,
                step_profile_phase=np.pi,
                length = 10,
                high_res=True
                )
random_road :flx.Road = flx.Road.make_random_road(length=10,
                                 smallest_wave_length=tyre_radius,
                                 frequency_scale=2,
                                 max_range=0.5)
fnt_road :flx.Road = flx.Road.make_road_from_file(Path.joinpath(repo_root_path, "data", "FinsAndThings_2d.asc"),
                                                  step_size = 0.01)
road = fnt_road
tyre = flx.ContinousTyre(initial_x=initial_x,
                          boundary_condition_file= matlab_file_path,
                          mass=unsprung_mass,
                          road=road,
                          free_radius=tyre_radius,
                          x_speed=forward_speed,
                          y_speed=0,
                          rigid_ring_nat_freq_hz=10,
                          rigid_ring_damping_ratio= 0.1,
                        )
q_car = flx.SprungMass(tyre_inst=tyre,
                       mass = sprung_mass,
                       speed_x=forward_speed,
                       speed_y=0,
                       spring_neutral_length=1,
                       natural_frequency_hz=1.5,
                       damping_ratio=0.4
                       )
data_logger : phsx.DataLogger = phsx.DataLogger()
data_logger.add_object(tyre)
data_logger.add_object(q_car)
data_logger.add_object(tyre.rigid_ring)
qmain_fig, Ax = plt.subplots(1 , 1)
tyre.find_new_contacts()
#for i in range(500): 
logged_data = []   
step = 0
current_ylim = (np.min(road.y) - 0.5 , np.max(road.y) + 0.5)
current_xlim = (road.x[0] , road.x[-1])
while tyre.states.position.x < road.length-2:
    step += 1
    plt.sca(Ax)
    st = time.time() # For timing the main operations
    '''
    main dynamics updates, should really be done in a function
    but because of the way the external forces are implemented, it's done 
    explicitly here
    '''
    q_car.update_states()
    tyre.update_states(-(q_car.spring_force + q_car.damper_force))
    # if len(tyre.contacts)>0:
    #     tyre.get_full_profile()
    # draw results
    if np.mod(step , draw_frequency) == 0:
        print(f'{1000*(time.time() - st):.1f} ms/t {q_car.states.velocity.y:0.3f}')
        # for ax in Ax:
        #     ax.cla()
        # Ax.cla()
        road.draw()
        tyre.draw()
        q_car.draw()
        plt.gca().set_aspect('equal')
        plt.sca(Ax)
        # for c in tyre.contacts:
        #     c.draw_pressure()
        plt.xlim(current_xlim)
        plt.ylim(current_ylim)
        plt.pause(0.001)
        current_ylim = plt.gca().get_ylim()
        current_xlim = plt.gca().get_xlim()

    if DEBUG:
        while not plt.waitforbuttonpress():
            pass
        current_ylim = plt.gca().get_ylim()
        current_xlim = plt.gca().get_xlim()
    Ax.cla()
data_logger.write_to_file(file_name=output_file_path)
print("SIMULTAION COMPLETE")


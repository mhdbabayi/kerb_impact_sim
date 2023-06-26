import flexring as flx
import physics_engine as phsx
import time
import os, sys
import numpy as np
from matplotlib import pyplot as plt
from euclid3 import Vector2
os.system("clear")
from scipy import io
# constants and initial values
DEBUG = sys.gettrace()
matlab_file_path = os.path.abspath("./kerb_impact_sim/data/beta_5.mat")
forward_speed_kph = 12.
forward_speed = forward_speed_kph/3.6
tyre_radius = 0.788/2
if DEBUG:
    draw_frequency = 1
    initial_x = 2
else:
    draw_frequency = 50
    initial_x = 1

sprung_mass = 700
unsprung_mass = 50
initial_y = tyre_radius - (sprung_mass + unsprung_mass)*10/(flx.ContinousTyre.lump_stiffness)
# defining sim objects, all moving objects inherit from rigid body
road = flx.Road(
                step_width=0.05,
                step_height=0.2,
                step_profile_phase=np.pi,
                high_res=True
                )
tyre = flx.ContinousTyre(initial_x=initial_x,
                          initial_y=initial_y+0.03,
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
while tyre.states.position.x < 4:
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
    if len(tyre.contacts)>0:
        tyre.get_full_profile()
    # draw results
    if np.mod(step , draw_frequency) == 0:
        print(f'{1000*(time.time() - st):.1f} ms/t {q_car.states.velocity.y:0.3f}')
        # for ax in Ax:
        #     ax.cla()
        # Ax.cla()
        tyre.draw()
        q_car.draw()
        plt.plot(road.x, road.y, color="brown")
        plt.gca().set_aspect('equal')
        plt.xlim(tyre.states.position.x+1.2*tyre.free_radius*\
                np.array((-1., 1.)))
        plt.ylim((-0.1, q_car.states.position.y + tyre.free_radius))
        plt.sca(Ax)
        # for c in tyre.contacts:
        #     c.draw_pressure()
        plt.pause(0.001)

    if DEBUG:
        while not plt.waitforbuttonpress():
            pass
    #[plt.sca(ax) for ax in Ax]
    Ax.cla()
for ax in Ax:
    ax.cla()

data_logger.write_to_file(file_name=os.path.abspath("./kerb_impact_sim/data/step_out.mat"))
print("SIMULTAION COMPLETE")


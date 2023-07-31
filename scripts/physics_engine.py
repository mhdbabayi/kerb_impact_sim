import numpy as np
from enum import Enum
from euclid3 import Vector2
from dataclasses import dataclass
from scipy import io
import matplotlib.pyplot as plt
import matplotlib.patches as patches
class ConstraintType(Enum):
    '''
    Name of constraint referes to the axis that's not updated 
    in the update method. The axis can still be explicitly changed 
    '''
    X = 1 << 0
    Y= 1 << 1
    THETA = 1 <<2
    

class DynamicObject:
    simParameters = {
    "universal_time" : 0,
    "time_step" : 0.001}
    graphic_objects = []
    def initialize(self):
        pass
    def iterate(self):
        pass
    @staticmethod
    def plot(x , y, *args, **kwargs):
        DynamicObject.graphic_objects += plt.plot(x , y )
    @staticmethod
    def add_patch(patch, *args, **kwargs):
        DynamicObject.graphic_objects.append(plt.gca().add_patch(patch, *args, **kwargs))
    def clear_plot():
        while len(DynamicObject.graphic_objects) > 0 :
            g = DynamicObject.graphic_objects.pop()
            plt.gca().get_children()[plt.gca().get_children().index(g)].remove()
            del g
        plt.draw()
class RigidBody(DynamicObject):
    def __init__(self,
                mass,
                initial_x,
                initial_y,
                initial_x_dot,
                initial_y_dot,
                constraint_type:str,
                name:str,
                initial_force_y = 0,
                initial_force_x = 0) -> None:
        self.states = RigidBody.State(
                                      position=Vector2(initial_x, initial_y),
                                      velocity=Vector2(initial_x_dot, initial_y_dot),
                                      acceleration= Vector2(0 , 0))
        self.forces = Vector2(initial_force_x , initial_force_y)
        self.mass = mass
        self.data_logger : DataLogger = None
        self.name = name
        '''
        the constraint is a string with three digits, each representing a bit
        from left to right, Theta, Y, X
        for example, '001' means, object is constrained in the X direction, and theta and Y are updated in the dynamics 
        '''
        self.constraint = constraint_type
    def update_states(self, external_forces :Vector2 =Vector2(0 , 0)):
        self.update_forces(external_forces)
        self.update_accelerations()
        self.states.velocity = self.states.velocity + self.states.acceleration * self.simParameters["time_step"]
        self.states.position = self.states.position + self.states.velocity * self.simParameters["time_step"]
        self.log_states()
    def update_accelerations(self):
        if not(int(self.constraint) & ConstraintType.X.value):
            self.states.acceleration.x = self.forces.x / self.mass
        if not(int(self.constraint) & ConstraintType.Y.value):
            self.states.acceleration.y = self.forces.y/self.mass
    def update_forces(self, external_forces:Vector2 = Vector2(0 , 0)):
        self.forces = Vector2(0, 0)
        self.forces = self.forces + external_forces
    def iterate(self):
        self.update_forces()
        self.update_accelerations()
        self.update_states()
    def log_states(self):
        if self.data_logger is not None:
            self.data_logger.data_dict[self.name]["position"]["x"].append(
                self.states.position[0])
            self.data_logger.data_dict[self.name]["position"]["y"].append(
                self.states.position[1])
            self.data_logger.data_dict[self.name]["velocity"]["x"].append(
                self.states.velocity[0])
            self.data_logger.data_dict[self.name]["velocity"]["y"].append(
                self.states.velocity[1])
    def log_force(self, force_value:Vector2):
        self.data_logger.data_dict[self.name]["force"]["x"].append(force_value[0])
        self.data_logger.data_dict[self.name]["force"]["y"].append(force_value[1])

    @dataclass
    class State:
        position = Vector2(0. , 0.)
        velocity = Vector2(0., 0.)
        acceleration = Vector2(0., 0.)
        def __init__(self, position, velocity, acceleration):
            self.position = position
            self.velocity = velocity
            self.acceleration = acceleration
            
class DataLogger(DynamicObject):
    def __init__(self) -> None:
        self.data_dict = {}
    def add_object(self, obj:RigidBody):
        if obj.name in self.data_dict:
            return
        obj.data_logger = self
        self.data_dict[obj.name] = {"position":{"x":[], "y":[]},
                                    "velocity": {"x":[], "y":[]},
                                    "force":{"x":[], "y":[]}}
    def write_to_file(self, file_name:str):
        for obj_name in self.data_dict.keys():
            for signal_name in self.data_dict[obj_name].keys():
                self.data_dict[obj_name][signal_name]["x"] = np.array(
                    self.data_dict[obj_name][signal_name]["x"])
                self.data_dict[obj_name][signal_name]["y"] = np.array(
                    self.data_dict[obj_name][signal_name]["y"])
        io.savemat(file_name, self.data_dict)
        


    



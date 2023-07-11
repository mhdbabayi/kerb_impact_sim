from typing import Any
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import physics_engine as phsx
from euclid3 import Vector2
from dataclasses import dataclass
import math_utils as ut
from scipy import interpolate
from scipy import io
import bisect
'''
sign convention: 
Tyre node zero is at the top and the nodes go counter-clockwise
Deflection of a node is positive towards the outside of the tyre(increase in radius)
'''



   
def fit_poly(P1, P2, P3):
    x1, y1, dy1 = (P1[0], P1[1], P1[2])
    x2, y2, dy2 = (P2[0], P2[1], P2[2])
    x3, y3, dy3 = (P3[0], P3[1], P3[2])

    A = np.array([
        [x1**5, x1**4, x1**3, x1**2, x1, 1],
        [5*x1**4, 4*x1**3, 3*x1**2, 2*x1, 1, 0],
        [x2**5, x2**4, x2**3, x2**2, x2, 1],
        [5*x2**4, 4*x2**3, 3*x2**2, 2*x2, 1, 0],
        [x3**5, x3**4, x3**3, x3**2, x3, 1],
        [5*x3**4, 4*x3**3, 3*x3**2, 2*x3, 1, 0]
    ])

    B = np.array([y1, dy1, y2, dy2, y3, dy3])
    coefficients = np.linalg.solve(A, B)
    
    return np.poly1d(coefficients)
def construct_piecewise_poly(start, end, peak):
    x1, y1 = (start[0],start[1])
    x2, y2 = (end[0], end[1])
    x3, y3 = (peak[0], peak[1])

    a1 = (y1-y3)/((x1-x3)**2)
    a2 = (y2 -y3)/(x2**2 - 2*x2*x3 + x3**2)
    b1 = -(2*x3*y1 - 2*x3*y3)/((x1-x3)**2)
    b2 = -(2*x3*y2 - 2*x3*y3)/(x2**2 - 2*x2*x3 + x3**2)
    c1 = (y3*x1**2 - 2*y3*x1*x3 + y1*x3**2)/((x1 - x3)**2)
    c2 = (y3*x2**2 - 2*y3*x2*x3 + y2*x3**2)/(x2**2 - 2*x2*x3 + x3**2)

    def piecewise_polynomial(x):
        if x <= x3:
            return a1*x**2 + b1*x + c1
        else:
            return a2*x**2 + b2*x + c2

    return piecewise_polynomial
def interpolate_boundary_condition(centre:Vector2,
                                   radius:float,
                                   point1:Vector2,
                                   point2:Vector2,
                                   dydx:float,
                                   ddydx:float,
                                   beta:float,
                                   direction:int):
    # given a chord or a semi chord on a circle, finds the 
    # point along the chord where separation happens
    dr_1 = (point1 - centre).magnitude() - radius
    dr_2 = (point2 - centre).magnitude() - radius
    dr_dtheta_1 = polar_derivative(point= point1 - centre,dy=dydx)
    dr_dtheta_2 = polar_derivative(point= point2- centre,dy= dydx)
    ddr_dtheta_1 = polar_second_derivative(point1 - centre,dydx,ddydx)
    ddr_dtheta_2 = polar_second_derivative(point2 - centre, dydx,ddydx)                                           
    # condition is true when lhs > rhs
    # we calculate lhs and rhs at point1 and point2 and linearly interpolate
    # to find the point where they cross
    """
    print(f'''r1: {dr_1}
    r2: {dr_2}\ndr1:{dr_dtheta_1}
    dr2: {dr_dtheta_1}
    ddr1: {ddr_dtheta_1}
    ddr2: {ddr_dtheta_2}''')
    """
    lhs_1 = 0.5*ddr_dtheta_1
    lhs_2 = 0.5*ddr_dtheta_2
    rhs_1 = -2*(beta**2)*(direction*dr_dtheta_1/beta + dr_1)
    rhs_2 = -2*(beta**2)*(direction*dr_dtheta_2/beta + dr_2)
    t = (rhs_1 - lhs_1)/((lhs_2-lhs_1)-(rhs_2 - rhs_1))# where rhs == lhs
    return point1 + t*(point2 - point1), t

class Road:
    def __init__(self, step_width, step_height,step_profile_phase = np.pi, length = 5, high_res=False) -> None:
        self.length = length
        self.x , self.y = self.create_profile(step_width, step_height , step_profile_phase, length)
        self.points = [Vector2(x , y) for x, y in zip(self.x , self.y)]
        self.high_res = high_res
        if high_res:
            self.points = self.over_sample(self.points)
            self.x = [p.x for p in self.points]
            self.y = [p.y for p in self.points]
        self.dydx = np.zeros(len(self.x))
        self.ddydx = np.zeros(len(self.x))
        for i in np.arange(start=1,stop=(len(self.x))-1):
            self.dydx[i] = (self.y[i+1] - self.y[i-1])/(self.x[i+1] - self.x[i-1])
            self.ddydx[i] = (self.y[i+1] + self.y[i-1] - 2*self.y[i])/\
                               ((self.x[i+1]-self.x[i])*(self.x[i] - self.x[i-1])) 
    @staticmethod
    def create_profile(step_width, step_height, step_profile_phase, length):
        x1 = np.array([0, length/2- step_width/2])
        y1 = np.array([0,0])
        x_step = length/2 + np.linspace(-step_width/2 , step_width/2,
                                        np.max((20,np.int32(step_width/0.01))))[1:]
        step_phase = np.linspace(0, step_profile_phase, len(x_step))
        y_step = step_height/2 - (step_height/2)*np.cos(step_phase)
        x2 = np.array([x_step[-1] + 0.01, length])
        y2 = np.array([y_step[-1], y_step[-1]])
        return np.hstack((x1 , x_step , x2)), np.hstack((y1, y_step , y2))
    @staticmethod
    def over_sample(points, distance_step = 0.001):
        cum_distance = [(points[idx+1] - points[idx]).magnitude()
                 for idx in range(len(points)-1)]
        cum_distance.insert(0 , 0)
        cum_distance = np.cumsum(cum_distance)
        x  = [p.x for p in points]
        y = [p.y for p in points]
        uniform_distance = np.arange(start=0, stop = cum_distance[-1], step=distance_step)
        x_interpolator = interpolate.interp1d(x = cum_distance, y = x,kind="linear")
        y_interpolator = interpolate.interp1d(x = cum_distance, y = y , kind="linear")
        return [Vector2(x_interpolator(d) , y_interpolator(d)) for d in uniform_distance]
class SmartRaod(Road):
    def __init__(self, step_width, step_height,step_profile_phase = np.pi, length = 5) -> None:
        # for now we initialize like before, with only a sine bump in the middle
        super().__init__(step_width, step_height,step_profile_phase, length)
        self.over_sampled_points = self.over_sample(self.points)
        self.node_list = []

    def initialize_nodes(self):
        self.node_list = [SmartRaod.Node(parent_road=self,
                                         position=self.points[idx],
                                         dydx=self.dydx[idx],
                                         ddydx=self.ddydx[idx]) for idx in range(len(self.points))]
        self.node_list[0].next = self.node_list[1]
        self.node_list[-1].prev = self.node_list[-2]
        for i in range(len(self.node_list)-2):
            self.node_list[i+1].next = self.node_list[i+2]
            self.node_list[i+1].prev = self.node_list[i]
        
    class Node():
        def __init__(self,
                     parent_road,
                     position:Vector2,
                     dydx:float,
                     ddydx:float,
                     next=None,
                     prev=None) -> None:
            self.parent_road = parent_road
            self.posistion:Vector2 = position
            self.dydx = dydx
            self.ddydx = ddydx
            self.next:SmartRaod.Node = next
            self.prev:SmartRaod.Node = prev

class ContinousTyre(phsx.RigidBody):
    beta = 5
    lump_stiffness = 500000.
    lump_damping = 500
    element_stiffness = 600000
    element_damping = 1000
    
    def __init__(self, initial_x, initial_y,road:Road,
                 boundary_condition_file:str,
                 free_radius = 1., node_res_deg = 1.,
                 x_speed = 0, y_speed = 0, mass = 50,
                 rigid_ring_nat_freq_hz = 100,
                 rigid_ring_damping_ratio = 0.5) -> None:
        super().__init__(mass = mass, initial_x=initial_x, initial_y = initial_y,name="tyre",
                       initial_x_dot = x_speed, initial_y_dot = y_speed,constraint_type='101'
                       )
        self.road = road
        self.free_radius = free_radius
        self.collisions = []
        self.contacts = []
        self.total_contact_force = Vector2(0 , 0)
        # nodes are only used for visualisation
        self.node_angles = np.deg2rad(np.linspace(0 , 360, 361)[0:-1])
        self.external_forces = 0
        self.rigid_ring =ContinousTyre.RigidRing(tyre = self,
                                                 natural_freq_hz=rigid_ring_nat_freq_hz,
                                                 damping_ratio=rigid_ring_damping_ratio,
                                                 )
        self.beam = ut.BeamTyre(beta = self.beta,
                                tyre_radius=self.free_radius,
                                boundary_theta_map_file=boundary_condition_file)
    def find_new_collisions(self, start_idx=0):
        road_idx = start_idx
        while self.road.points[road_idx].x < self.states.position.x - self.free_radius:
                road_idx += 1
        while self.road.points[road_idx].x < self.states.position.x + self.free_radius:
            T= ut.circle_line_intersection(self.road.points[road_idx],
                                              self.road.points[road_idx+1],
                                              self.states.position,
                                              self.free_radius) 
            if T is not None:
                if T[0] is not None:
                    self.collisions.append(ContinousTyre.Collision(start=T[0],
                                                    end=None,
                                                    start_road_idx = road_idx,
                                                    end_road_idx = road_idx))
                    # if the line crosses in only 1 point
                    while self.collisions[-1].fore_point is None:
                        road_idx +=1
                        T = ut.circle_line_intersection(self.road.points[road_idx],
                                                        self.road.points[road_idx+1],
                                                        self.states.position,
                                                        self.free_radius)
                        if T is not None:
                            self.collisions[-1].fore_point = T[1]
                            self.collisions[-1].end_road_idx = road_idx            
            road_idx +=1           
    def find_new_contacts(self, start_idx = 0):
        self.find_new_collisions(start_idx)
        for c in self.collisions:
            #self.contacts.append(ContinousTyre.Contact(self, c))
            bisect.insort(self.contacts, ContinousTyre.Contact(self, c))
            self.collisions.remove(c)
    def draw(self):
        circle_obj = patches.Circle(self.states.position, self.free_radius, linewidth=1,fill=False)
        plt.gca().add_patch(circle_obj)
        for c in self.contacts:
            c.draw()
        self.rigid_ring.draw()
    def update_forces(self, external_forces):
        super().update_forces(external_forces)
        self.total_contact_force = Vector2(0 , 0)
        for c in self.contacts:
            self.total_contact_force = self.total_contact_force + c.get_forces_centre_point()
        self.forces += self.total_contact_force
        self.log_force(self.total_contact_force)
    def update_states(self, external_forces):
        self.external_forces = external_forces
        super().update_states(external_forces)
        self.update_contacts()
        self.rigid_ring.update_states()
    def update_contacts(self):
        if len(self.contacts) == 0:
            start_idx = 0
        else:
            start_idx=self.contacts[-1].collision.end_road_idx+5
        self.find_new_contacts(start_idx)
        self.contacts = [c for c in self.contacts if c.update() is not None]    
    def get_full_profile(self):
        theta = np.arange(0, np.deg2rad(360), np.deg2rad(1))
        i = 0
        while i < len(self.contacts)-1:
            print(np.arange(np.rad2deg(self.contacts[i].aft_theta_abs()[0]),
                            np.rad2deg(self.contacts[i+1].fore_theta_abs()[0]),
                            1))
            i+=1
        print(np.arange(np.rad2deg(self.contacts[i].fore_theta_abs()[0]),
                            np.rad2deg(self.contacts[-1].aft_theta_abs()[0])+360,
                            1))
        print("**************")



    @dataclass
    class Collision:
        aft_point:Vector2
        fore_point:Vector2
        road_idx:int
        def __init__(self, start, end, start_road_idx, end_road_idx):
            self.aft_point = start
            self.fore_point = end
            self.start_road_idx = start_road_idx
            self.end_road_idx = end_road_idx
        def update(self, centre_road_idx, tyre_centre, tyre_radius,road_inst):
            self.end_road_idx = centre_road_idx
            self.start_road_idx = centre_road_idx
            while ((tyre_centre - road_inst.points[self.end_road_idx]).magnitude() < \
                tyre_radius):
                self.end_road_idx += 1
            self.end = road_inst.points[self.end_road_idx]
            while ((tyre_centre - road_inst.points[self.start_road_idx]).magnitude() < \
                tyre_radius):
                self.start_road_idx -= 1
            self.aft_point = road_inst.points[self.start_road_idx]          
    class Contact:
        def __init__(self,
                     tyre,
                     collision):
            self.tyre: ContinousTyre = tyre
            self.collision: ContinousTyre.Collision = collision
            self.centre_point_idx = np.argmin((p - self.tyre_states.position).magnitude() for p in\
                                               self.tyre.road.points[collision.start_road_idx:collision.end_road_idx])+\
                                                  collision.start_road_idx
            self.centre_point = lambda : self.tyre.road.points[self.centre_point_idx]
            self.centre_point_angle = lambda : np.arctan2(self.centre_point().y - self.tyre.states.position.y,
                                                 self.centre_point().x - self.tyre.states.position.x)
            self.centre_point_deformation = lambda : self.tyre.free_radius - (self.tyre.states.position - self.centre_point()).magnitude()
            self.normal_vector = lambda : self.centre_point() - self.tyre.states.position

            self.fore_theta_profile = None
            self.aft_theta_profile = None
            self.fore_theta_abs = lambda : self.fore_theta_profile + self.centre_point_angle()
            self.aft_theta_abs = lambda : -self.aft_theta_profile + self.centre_point_angle()
            self.fore_deformation_profile = None
            self.aft_deformation_profile = None

            self.fore_circle_centre = None
            self.fore_circle_radius = None
            self.aft_circle_centre = None
            self.aft_circle_radius = None

            self.fit_theta = None
            self.fit_deformation = None
            self.fit_theta_old = None
            self.fit_deformation_old = None

            self.whole_theta_profile = None
            self.whole_deformation_profile = None
            self.prev_whole_deformation_profile = None
            self.prev_whole_theta_profile = None


            #hack
            self.normal_force = 0
            self.centre_migration_theta = None # movement of contact centre relative to previous step
            self.prev_centre_point_deformation = 0

            self.update()
        def __lt__(self, other):
            if not isinstance(other, ContinousTyre.Contact):
                return NotImplemented
            return self.centre_point_angle() < other.centre_point_angle()
        def draw(self):
            plt.plot(self.centre_point().x , self.centre_point().y , "o")
            self.draw_terrain_circles()
            self.draw_envelop()
        def draw_pressure(self):
            # plt.plot(np.rad2deg(self.prev_whole_theta_profile),
            #          self.prev_whole_deformation_profile*1000 ,"r--")
            # plt.plot(np.rad2deg(self.whole_deformation_profile),
            #          self.whole_deformation_profile*1000 , 'b--')
            # plt.xlabel("theta deg") 
            # plt.ylabel("sidewall element deformation mm")
            plt.plot(self.tyre.data_logger.data_dict["tyre"]["position"]["x"],
            self.tyre.data_logger.data_dict["tyre"]["force"]["y"])          
        def draw_terrain_circles(self):
            fore_circle = patches.Wedge(center=self.fore_circle_centre,
                                        r=self.fore_circle_radius,
                                        theta1= np.rad2deg(self.centre_point_angle()),
                                        theta2= np.rad2deg(self.centre_point_angle()) + 180,
                                        fill = None,
                                        width=0.0001,
                                        lw = 2,
                                        color = "magenta",
                                        linestyle = "--"
                                         )
            aft_circle = patches.Wedge(center = self.aft_circle_centre,
                                       r= self.aft_circle_radius,
                                       theta1 = np.rad2deg(self.centre_point_angle()) - 180,
                                       theta2 = np.rad2deg(self.centre_point_angle()),
                                       fill= None,
                                       width = 0.0001,
                                       lw =2,
                                       color = "magenta",
                                       linestyle = "--")
            plt.gca().add_patch(fore_circle)
            plt.gca().add_patch(aft_circle)
        def draw_envelop(self):
            x = [self.tyre.states.position.x +\
                  (self.tyre.free_radius - w)*np.cos(t) for
                    w,t in zip(self.whole_deformation_profile, self.whole_theta_profile)]
            y = [self.tyre.states.position.y +\
                  (self.tyre.free_radius - w)*np.sin(t) for
                    w,t in zip(self.whole_deformation_profile, self.whole_theta_profile)]
            plt.plot(x , y , "green")
        def set_equivalent_circles(self):
            # TODO sign of tangent input acts strange
            fore_curvature = ut.get_circle_tangent_2points(tangent=-self.normal_vector().cross(),
                                                           p0= self.centre_point(),
                                                           p1= self.tyre.road.points[self.centre_point_idx + 15])
            aft_curvature = ut.get_circle_tangent_2points(tangent=self.normal_vector().cross(),
                                                          p0 = self.centre_point(),
                                                          p1= self.tyre.road.points[self.centre_point_idx - 15])
            # check for zero curvature
            if np.abs(fore_curvature) < 0.1:
                self.fore_circle_radius = 10
            else: 
                self.fore_circle_radius = 1/fore_curvature
            if np.abs(aft_curvature) < 0.1:
                self.aft_circle_radius = 10
            else:
                self.aft_circle_radius = 1/aft_curvature
            self.fore_circle_centre = self.centre_point() +\
                self.fore_circle_radius*self.normal_vector().normalized()
            self.aft_circle_centre= self.centre_point() +\
                self.aft_circle_radius*self.normal_vector().normalized()
        def get_forces_centre_point(self):
            spring_force = self.tyre.lump_stiffness * self.centre_point_deformation()
            damping_force = self.tyre.lump_damping *\
                (self.centre_point_deformation() -\
                  self.prev_centre_point_deformation)/self.tyre.simParameters["time_step"] 
            self.normal_force = spring_force + damping_force
            #print(f'\tspring force: {spring_force:.1f}\n')
            return -self.normal_force*self.normal_vector().normalized()
        def set_deformation(self):
            penetration = self.centre_point_deformation()
            self.fore_theta_profile, self.fore_deformation_profile =\
                self.tyre.beam(penetration, self.fore_circle_radius)
            self.aft_theta_profile , self.aft_deformation_profile =\
                self.tyre.beam(penetration, self.aft_circle_radius)                             
        def get_stacked_profile(self):
            # returns the full profile of the contact, as two vectors, theta and w
            # each 121 elements, covering 120 degrees
            angle_step = self.tyre.beam.theta_resolution 
            contact_section_theta = np.arange(start = self.aft_theta_abs()[0] + angle_step,
                                              stop= self.fore_theta_abs()[0],
                                              step = angle_step) - self.centre_point_angle()
            contact_section_deformation = ut.fit_quadratic(left_point = (contact_section_theta[0], self.aft_deformation_profile[0]),
                                                           right_point= (contact_section_theta[-1], self.fore_deformation_profile[0]),
                                                           y0 = self.centre_point_deformation())(contact_section_theta)
            stacked_theta = np.hstack((np.flip(self.aft_theta_abs()), contact_section_theta + self.centre_point_angle(), self.fore_theta_abs()))
            stacked_deformation = np.hstack((np.flip(self.aft_deformation_profile),
                                         contact_section_deformation,
                                         self.fore_deformation_profile))
            uniform_theta = np.linspace(np.deg2rad(-60) , np.deg2rad(60) , 121) + self.centre_point_angle()
            uniform_deformation = interpolate.interp1d(stacked_theta , stacked_deformation)(uniform_theta)
            return  uniform_theta, uniform_deformation
        def update(self):
            self.update_centre_point_idx() 
            # check if contact still exists
            self.set_equivalent_circles()
            if self.centre_point_deformation() < 0:
                return None
            self.collision.update(self.centre_point_idx,
                                  self.tyre.states.position,
                                  self.tyre.free_radius,
                                  self.tyre.road)
            self.update_deformation()
            return 1
        def update_centre_point_idx(self):
            #which way to go on the road to find the new centre point idx
            road_tangent_vector = self.tyre.road.points[self.centre_point_idx+1] -\
                self.tyre.road.points[self.centre_point_idx-1]
            idx_increment = np.int32(np.sign(
                    road_tangent_vector.dot(self.tyre.states.velocity)))
            prev_centre_idx = self.centre_point_idx
            self.prev_centre_point_deformation = self.centre_point_deformation()
            while ((self.tyre.road.points[self.centre_point_idx] - self.tyre.states.position).magnitude() >\
                (self.tyre.road.points[self.centre_point_idx+idx_increment] - self.tyre.states.position).magnitude()):
                self.centre_point_idx += idx_increment
            rolling_radius = (self.centre_point() - self.tyre.states.position).magnitude()
            self.centre_migration_theta = (self.tyre.road.points[prev_centre_idx] - self.centre_point()).magnitude()/rolling_radius
        def update_deformation(self):
            # TODO fix for first iteration
            if self.whole_theta_profile is not None:
                self.prev_whole_deformation_profile = self.whole_deformation_profile
                self.prev_whole_theta_profile = self.whole_theta_profile - self.centre_migration_theta
            self.set_deformation()
            self.whole_theta_profile , self.whole_deformation_profile = self.get_stacked_profile()
            if self.prev_whole_theta_profile is None:
                self.prev_whole_deformation_profile = self.whole_deformation_profile
                self.prev_whole_theta_profile = self.whole_theta_profile - self.centre_migration_theta
    class RigidRing(phsx.RigidBody):
        def __init__(self,tyre,
                     natural_freq_hz= 100,
                      damping_ratio = 0.5):
            super().__init__(mass= 10,
                            initial_x= tyre.states.position.x,
                            initial_y=tyre.states.position.y,
                            initial_x_dot=tyre.states.velocity.x,
                            initial_y_dot=tyre.states.velocity.y,
                            constraint_type= '100', name="rigid_ring")
            self.tyre:ContinousTyre = tyre
            self.stiffness = 2*np.pi * natural_freq_hz**2*self.mass
            self.damping = 2*2*np.pi*natural_freq_hz* self.mass * damping_ratio
            self.displacement_vector = lambda :\
                self.states.position - self.tyre.states.position
            self.velocity_vector = lambda :\
                self.states.velocity - self.tyre.states.velocity
            self.prev_displacement_vector = self.displacement_vector()
        def update_forces(self, external_forces:Vector2 = Vector2(0 , 0)):
            self.forces = -self.stiffness * self.displacement_vector() -\
                            self.damping * self.velocity_vector()
            self.log_force(self.forces)
        def draw(self):
            patch = patches.Circle(xy =(self.states.position),
                                   radius =0.8*self.tyre.free_radius,
                                   fill = False,
                                   lw= 2, 
                                   color = "red")
            plt.gca().add_patch(patch)
class SprungMass(phsx.RigidBody):
    def __init__(self,
                tyre_inst:ContinousTyre,
                mass,
                speed_x=0,
                speed_y=0,
                spring_neutral_length = 0.3,
                natural_frequency_hz = 1.5,
                damping_ratio = 0.5):
        self.tyre_inst = tyre_inst
        super().__init__(mass=mass, initial_x=tyre_inst.states.position.x,
                        initial_y = tyre_inst.states.position.y + spring_neutral_length,
                            initial_x_dot = self.tyre_inst.states.velocity.x,
                              initial_y_dot = speed_y,
                            constraint_type = '001', name="sprung_mass")
        self.natural_freq_rad = 360*np.deg2rad(natural_frequency_hz)
        self.damping_ratio = damping_ratio
        self.spring_neutral_length = spring_neutral_length
        self.spring_stiffness = (self.natural_freq_rad)**2 * self.mass
        self.damping_coefficient = 2*self.natural_freq_rad * self.mass * damping_ratio
        self.spring_preload = self.mass * 9.8
        self.spring_force = Vector2(0 ,self.spring_preload) 
        self.damper_force = Vector2(0 , 0)
    def update_forces(self, external_forces:Vector2 = Vector2(0 ,0)):
        self.spring_force = Vector2(0 ,(self.spring_neutral_length - 
        (self.states.position.y - self.tyre_inst.states.position.y)) + self.spring_preload)
        self.damper_force = Vector2(0,(self.tyre_inst.states.velocity.y  - self.states.velocity.y) * self.damping_coefficient)
        gravity_force = Vector2(0 , -self.mass*9.8)
        self.forces =  self.spring_force + self.damper_force + gravity_force       
        self.log_force(self.spring_force + self.damper_force)
    def draw(self):
        width = self.tyre_inst.free_radius * 2
        rect = patches.Rectangle((self.states.position.x-width/2, self.states.position.y-width/2),
                                 width=width, height=width/2)
        plt.gca().add_patch(rect)
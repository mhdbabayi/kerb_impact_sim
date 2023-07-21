import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import physics_engine as phsx
from euclid3 import Vector2
from dataclasses import dataclass
import math_utils as ut
from scipy import interpolate
from scipy import io
import bisect
from pathlib import Path
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
    filter_window_size_global = 20 
    def __init__(self, length = 5, high_res=False) -> None:
        self.x = None
        self.y = None
        self.dydx = None
        self.ddydx = None
        self.points:list[Vector2] = None
        self.length = length
        self.high_res = high_res
        self.random_profile_seed = None
        self.start_node:Road.Node = None
        self.start_section : Road.Section = None
    def setup_points(self):
        self.points = [Vector2(x , y) for x, y in zip(self.x , self.y)]
        if self.high_res:
            self.points = self.over_sample(self.points)
            self.x = np.array([p.x.item() for p in self.points])
            self.y = np.array([p.y.item() for p in self.points])
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
    @staticmethod
    def make_simple_road(step_width, step_height, step_profile_phase= np.pi, length = 5, high_res = False):
        road = Road(length=length,high_res=high_res)
        road.x , road.y = Road.create_profile(step_width, step_height , step_profile_phase, length)
        road.setup_points()
        road.make_smart()
        return road
    @staticmethod
    def make_random_road(length,
                         smallest_wave_length,
                         frequency_scale = 1,
                         seed= None,
                         max_range = None):
        road = Road(length=length, high_res=True)
        road.x = np.arange(0 , length, 0.01)
        road.y = np.zeros(len(road.x))
        road.random_profile_seed = seed
        frequencies = 2*np.pi/(np.linspace(smallest_wave_length, length, 10))
        rng = np.random.default_rng(seed=seed)
        weights = rng.random((len(frequencies)))/(np.log(frequency_scale*frequencies))
        phases = rng.random((len(frequencies)))*2*np.pi
        for idx in range(len(frequencies)):
            road.y = road.y + weights[idx]*np.sin(road.x *frequencies[idx] + phases[idx])
        if max_range is not None:
            road.y = road.y *max_range/ (np.max(road.y) - np.min(road.y))
        road.setup_points()
        road.make_smart()
        return road
    @staticmethod
    def make_road_from_file(file_name, step_size=0.001):
        with open(file_name) as f:
            point_list = np.array([[float(f) for f in line.split(",")] for line in f])
        x , y = Road.under_sample(x = point_list[: , 0] , z=point_list[: , 2] , x_step=step_size)
        road = Road(length=x[-1], high_res=False)
        road.x = x
        road.y = y
        road.setup_points()
        road.make_smart()
        return road    
    @staticmethod
    def under_sample(x: np.array, z:np.array, x_step):
        # under sample x and z by picking a mean z in every x bin with width x_step
        xbins = np.digitize(x , np.arange(start=x[0], stop=x[-1], step=x_step))
        out_z = []
        out_x = []
        for idx,x_sample in enumerate(np.arange(start=x[0], stop=x[-1], step = x_step)):
            out_z.append(np.mean(z[xbins == idx+1]))
            out_x.append(x_sample)
        return np.array(out_x)[~np.isnan(out_z)] , np.array(out_z)[~np.isnan(out_z)]
    def make_smart(self):
        self.start_node = Road.Node(idx=0, road = self)
        end_node = self.start_node
        while end_node.make_next() is not None:
            end_node = end_node.next
        self.start_section = Road.Section(road = self, number=0, start_node=self.start_node)
        end_section = self.start_section
        while end_section.make_next() is not None:
            end_section = end_section.next
        s = self.start_section
        while s is not None:
            if s.end_node is s.start_node:
                s.remove()
            s = s.next
        
    def draw(self):
        plt.plot(self.x , self.y , color="brown")
        s = self.start_section
        while s is not None:
            s.draw()
            s = s.next
    class Node():
        def __init__(self,
                     idx,
                     road,
                     prev= None) -> None:
            self.road:Road = road
            self.section:Road.Section = None
            self.position:Vector2 = road.points[idx]
            self.prev:Road.Node = prev
            self.next:Road.Node = None
            self.idx = idx
            if idx > self.road.filter_window_size_global and idx < (len(self.road.points) - self.road.filter_window_size_global):
                self.tangent = (self.road.points[idx+ self.road.filter_window_size_global] -\
                                self.road.points[idx - self.road.filter_window_size_global]).normalized()
                self.curvature = ut.get_curvature_3point(prev = self.road.points[idx - self.road.filter_window_size_global],
                                                        target= self.position,
                                                        next = self.road.points[idx + self.road.filter_window_size_global])
            else:
                self.tangent = Vector2(0 , 0)
                self.curvature = 0
        def make_next(self):
            if self.next is not None:
                return self.next
            if self.idx < len(self.road.points) -1:
                self.next = Road.Node(idx = self.idx+1, road = self.road, prev=self)
                return self.next
            return None
            
    class Section():
        def __init__(self,road,number, start_node, prev= None) -> None:
            self.road:Road = road
            self.number = number
            self.prev:Road.Section = prev
            self.next:Road.Section = None
            self.start_node:Road.Node = start_node
            self.max_curvature_mag = np.abs(start_node.curvature)
            self.end_node:Road.Node = start_node 
            self.peak_node:Road.Node = start_node
            self.curvature_sign = np.sign(start_node.curvature)
            while np.sign(self.end_node.next.curvature) == np.sign(self.end_node.curvature):
                self.end_node.section = self
                if np.abs(self.end_node.curvature) > self.max_curvature_mag:
                    self.peak_node = self.end_node
                    self.max_curvature_mag = np.abs(self.end_node.curvature)
                self.end_node = self.end_node.next
                if self.end_node.next is None:
                    break
        def make_next(self):
            if self.next is not None:
                return self.next    
            elif self.end_node.next is not None:
                self.next = Road.Section(road = self.road, number=self.number+1, start_node=self.end_node.next,
                                            prev=self)
                return self.next
            else:
                return None
        def draw(self):
            if self.curvature_sign > 0:
                color = "red"
            else:
                color = "black"
            points = np.array(self.road.points[self.start_node.idx:self.end_node.idx])
            if len(points>0):
                plt.plot(points[: , 0] , points[: , 1] , color= color)
                plt.plot(self.peak_node.position.x , self.peak_node.position.y , "m*")
        def remove(self):
            if self.next is not None:
                self.next.start_node = self.end_node
                self.next.number = self.number
                self.next.prev = self.prev
            if self.prev is not None:
                self.prev.end_node = self.start_node
                self.prev.next = self.next
        def check_contact_initial(self, tyre):
            # returns closest node in segment and the next segment if it's still in the bounding box
            if self.curvature_sign >0:
                return None, self.next
            min_node:Road.Node = None
            current_node = self.start_node
            min_distance_sqr = tyre.free_radius**2
            while current_node is not self.end_node:
                if (current_distance_sqr:=(current_node.position - tyre.states.position).magnitude_squared())\
                      < min_distance_sqr:
                    min_distance_sqr = current_distance_sqr
                    min_node = current_node
                current_node = current_node.next
            if  self.end_node.position.x > (tyre.states.position.x + tyre.free_radius):
                return min_node, None
            return min_node, self.next

class ContinousTyre(phsx.RigidBody):
    beta = 5
    lump_stiffness = 500000.
    lump_damping = 500
    element_stiffness = 1000000
    element_damping = 1000
    
    def __init__(self, initial_x,road:Road,
                 boundary_condition_file:str,
                 free_radius = 1., node_res_deg = 1.,
                 x_speed = 0, y_speed = 0, mass = 50,
                 rigid_ring_nat_freq_hz = 100,
                 rigid_ring_damping_ratio = 0.5) -> None:
        
        initial_y = road.y[np.where(road.x > initial_x)[0][0]] + free_radius *0.95
        super().__init__(mass = mass, initial_x=initial_x, initial_y = initial_y,name="tyre",
                       initial_x_dot = x_speed, initial_y_dot = y_speed,constraint_type='101'
                       )
        self.road = road
        self.free_radius = free_radius
        self.collisions:list[ContinousTyre.Collision] = []
        self.contacts:list[ContinousTyre.Contact] = []
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
        t0 = time.time()
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
        print(f"\t{(time.time()-t0)*1000:.1f} to find collisions")          
    def find_new_contacts(self, start_idx = 0):
        self.find_new_collisions(start_idx)
        for c in self.collisions:
            #self.contacts.append(ContinousTyre.Contact(self, c))
            bisect.insort(self.contacts, ContinousTyre.Contact(self, c))
            self.collisions.remove(c)
    def initialize_contact(self):
        section:Road.Section = self.road.start_section
        while section.end_node.position.x < (self.states.position.x - self.free_radius):
            section = section.next 
        contact_centres = []
        while section is not None:
            centre_node, section = section.check_contact_initial(tyre = self)
            if centre_node is not None:
                contact_centres.append(centre_node)
        return contact_centres
        
    def draw(self):
        circle_obj = patches.Circle(self.states.position, self.free_radius, linewidth=1,fill=False)
        phsx.DynamicObject.add_patch(circle_obj)
        for c in self.contacts:
            c.draw()
        self.rigid_ring.draw()
    def update_forces(self, external_forces):
        super().update_forces(external_forces)
        self.total_contact_force = Vector2(0 , 0)
        for c in self.contacts:
            #self.total_contact_force = self.total_contact_force + c.get_forces_centre_point()
            self.total_contact_force = self.total_contact_force + c.get_forces_integral()
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
            start_idx=self.contacts[-1].centre_point_idx + 10
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
        print(np.arange(np.rad2deg(self.contacts[i].fore_theta_abs_f()[0]),
                            np.rad2deg(self.contacts[-1].aft_theta_abs_f()[0])+360,
                            1))
        print("**************")
    @dataclass
    class Collision:
        aft_point:Vector2
        fore_point:Vector2
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
            self.centre_point_f = lambda : self.tyre.road.points[self.centre_point_idx]
            self.centre_point_angle_f = lambda : np.arctan2(self.centre_point_f().y - self.tyre.states.position.y,
                                                 self.centre_point_f().x - self.tyre.states.position.x)
            self.centre_point_deformation_f = lambda : self.tyre.free_radius - (self.tyre.states.position - self.centre_point_f()).magnitude()
            self.normal_vector_f = lambda : self.centre_point_f() - self.tyre.states.position

            self.fore_theta_profile = None # starting from fore separation_theta (relative to centre) (always positvie)
            self.aft_theta_profile = None # starting from aft separation theta and going up (always positive)
            self.fore_theta_abs_f = lambda : self.fore_theta_profile + self.centre_point_angle_f() # starting from fore separation theta and going up
            self.aft_theta_abs_f = lambda : -self.aft_theta_profile + self.centre_point_angle_f() # starting from aft separation theta and going DOWN
            self.fore_deformation_profile = None # element wise match with both profile and abs
            self.aft_deformation_profile = None # element wise match with both profile and abs

            self.fore_circle_centre : Vector2 = None # global coordinates
            self.fore_circle_radius = None
            self.aft_circle_centre : Vector2 = None # global coordinates
            self.aft_circle_radius = None

            self.fit_theta = None
            self.fit_deformation = None
            self.fit_theta_old = None
            self.fit_deformation_old = None

            self.whole_theta_profile = None # constant length and absolute
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
            return self.centre_point_angle_f() < other.centre_point_angle_f()
        def draw(self):
            graphic_objects = []
            graphic_objects.append(
                phsx.DynamicObject.plot(self.centre_point_f().x , self.centre_point_f().y , marker ="o")
            )
            self.draw_terrain_circles()
            self.draw_envelop()
        def draw_pressure(self):
            # phsx.DynamicObject.plot(np.rad2deg(self.prev_whole_theta_profile),
            #          self.prev_whole_deformation_profile*1000 ,"r--")
            # phsx.DynamicObject.plot(np.rad2deg(self.whole_deformation_profile),
            #          self.whole_deformation_profile*1000 , 'b--')
            # plt.xlabel("theta deg") 
            # plt.ylabel("sidewall element deformation mm")
            phsx.DynamicObject.plot(self.tyre.data_logger.data_dict["tyre"]["position"]["x"],
            self.tyre.data_logger.data_dict["tyre"]["force"]["y"])          
        def draw_terrain_circles(self):
            graphic_objects = []
            fore_circle = patches.Wedge(center=self.fore_circle_centre,
                                        r=self.fore_circle_radius,
                                        theta1= np.rad2deg(self.centre_point_angle_f()),
                                        theta2= np.rad2deg(self.centre_point_angle_f()) + 180,
                                        fill = None,
                                        width=0.0001,
                                        lw = 2,
                                        color = "magenta",
                                        linestyle = "--"
                                         )
            aft_circle = patches.Wedge(center = self.aft_circle_centre,
                                       r= self.aft_circle_radius,
                                       theta1 = np.rad2deg(self.centre_point_angle_f()) - 180,
                                       theta2 = np.rad2deg(self.centre_point_angle_f()),
                                       fill= None,
                                       width = 0.0001,
                                       lw =2,
                                       color = "magenta",
                                       linestyle = "--")
            graphic_objects.append(
                phsx.DynamicObject.add_patch(fore_circle)
            )
            graphic_objects.append(
                phsx.DynamicObject.add_patch(aft_circle)
            )
            return graphic_objects
        def draw_envelop(self):
            graphic_objects = []
            x = [self.tyre.states.position.x +\
                  (self.tyre.free_radius - w)*np.cos(t) for
                    w,t in zip(self.whole_deformation_profile, self.whole_theta_profile)]
            y = [self.tyre.states.position.y +\
                  (self.tyre.free_radius - w)*np.sin(t) for
                    w,t in zip(self.whole_deformation_profile, self.whole_theta_profile)]
            phsx.DynamicObject.plot(x , y , color = "green")
            return graphic_objects
        def set_equivalent_circles(self):
            # TODO sign of tangent input acts strange
            fore_curvature = ut.get_circle_tangent_2points(tangent=-self.normal_vector_f().cross(),
                                                           p0= self.centre_point_f(),
                                                           p1= self.tyre.road.points[self.centre_point_idx + 15])
            aft_curvature = ut.get_circle_tangent_2points(tangent=self.normal_vector_f().cross(),
                                                          p0 = self.centre_point_f(),
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
            self.fore_circle_centre = self.centre_point_f() +\
                self.fore_circle_radius*self.normal_vector_f().normalized()
            self.aft_circle_centre= self.centre_point_f() +\
                self.aft_circle_radius*self.normal_vector_f().normalized()
        def get_forces_centre_point(self):
            spring_force = self.tyre.lump_stiffness * self.centre_point_deformation_f()
            damping_force = self.tyre.lump_damping *\
                (self.centre_point_deformation_f() -\
                  self.prev_centre_point_deformation)/self.tyre.simParameters["time_step"] 
            self.normal_force = spring_force + damping_force
            #print(f'\tspring force: {spring_force:.1f}\n')
            return -self.normal_force*self.normal_vector_f().normalized()
        def get_forces_integral(self,
                                distance_to_previous_contact= None,
                                distance_to_next_contact=None) -> Vector2:
            contact_patch_forces = self.get_contact_patch_forces()
            fore_section_forces = self.tyre.element_stiffness*self.tyre.beam.get_normal_force_integral(
                w0 = self.fore_deformation_profile[0],
                dw0 = (self.fore_deformation_profile[1] - self.fore_deformation_profile[0])/\
                 (self.fore_theta_profile[1] - self.fore_theta_profile[0]), theta_end=distance_to_next_contact)
            aft_section_forces = self.tyre.element_stiffness*self.tyre.beam.get_normal_force_integral(
                w0 = self.fore_deformation_profile[0],
                dw0 = (self.fore_deformation_profile[1] - self.fore_deformation_profile[0])/\
                 (self.fore_theta_profile[1] - self.fore_theta_profile[0]), theta_end=distance_to_previous_contact)
            return (contact_patch_forces + fore_section_forces + aft_section_forces)*\
                -self.normal_vector_f().normalized()        
        def set_deformation(self):
            penetration = self.centre_point_deformation_f()
            self.fore_theta_profile, self.fore_deformation_profile =\
                self.tyre.beam(penetration, self.fore_circle_radius)
            self.aft_theta_profile , self.aft_deformation_profile =\
                self.tyre.beam(penetration, self.aft_circle_radius)                             
        def get_contact_patch_forces(self):
            # estimate contact section with two quadratic curves
            # of the form d0 - A(theta^2) wehre d0 in centre point deformation 
            # and deformation at separation is thet same as profile
            fore_A = (self.centre_point_deformation_f() - self.fore_deformation_profile[0])/\
                        self.fore_theta_profile[0]**2
            aft_A = (self.centre_point_deformation_f() - self.aft_deformation_profile[0])/\
                        self.aft_theta_profile[0]**2
            contact_spring_force = (self.centre_point_deformation_f() * (
                self.fore_theta_profile[0] + self.aft_theta_profile[0]) - (
                fore_A/3 * self.fore_theta_profile[0]**3 + 
                aft_A/3 * self.aft_theta_profile[0]**3))* self.tyre.element_stiffness
            return contact_spring_force
        def get_stacked_profile(self):
            # returns the full profile of the contact, as two vectors, theta and w
            # each 121 elements, covering 120 degrees
            angle_step = self.tyre.beam.theta_resolution 
            contact_section_theta = np.arange(start = self.aft_theta_abs_f()[0] + angle_step,
                                              stop= self.fore_theta_abs_f()[0],
                                              step = angle_step) - self.centre_point_angle_f()
            contact_section_deformation = ut.fit_quadratic(left_point = (contact_section_theta[0], self.aft_deformation_profile[0]),
                                                           right_point= (contact_section_theta[-1], self.fore_deformation_profile[0]),
                                                           y0 = self.centre_point_deformation_f())(contact_section_theta)
            stacked_theta = np.hstack((np.flip(self.aft_theta_abs_f()), contact_section_theta + self.centre_point_angle_f(), self.fore_theta_abs_f()))
            stacked_deformation = np.hstack((np.flip(self.aft_deformation_profile),
                                         contact_section_deformation,
                                         self.fore_deformation_profile))
            uniform_theta = np.linspace(np.deg2rad(-60) , np.deg2rad(60) , 121) + self.centre_point_angle_f()
            uniform_deformation = interpolate.interp1d(stacked_theta , stacked_deformation)(uniform_theta)
            return  uniform_theta, uniform_deformation
        def update(self):
            self.update_centre_point_idx() 
            # check if contact still exists
            self.set_equivalent_circles()
            if self.centre_point_deformation_f() < 0:
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
            self.prev_centre_point_deformation = self.centre_point_deformation_f()
            while ((self.tyre.road.points[self.centre_point_idx] - self.tyre.states.position).magnitude() >\
                (self.tyre.road.points[self.centre_point_idx+idx_increment] - self.tyre.states.position).magnitude()):
                self.centre_point_idx += idx_increment
            rolling_radius = (self.centre_point_f() - self.tyre.states.position).magnitude()
            self.centre_migration_theta = (self.tyre.road.points[prev_centre_idx] - self.centre_point_f()).magnitude()/rolling_radius
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
    class MultiContact:
        def __init__(self,
                     tyre, 
                     node_list:list[Road.Node],
                     ) -> None:
            self.tyre = tyre
            self.node_list: list[Road.Node] = sorted(node_list,
                                                      key=lambda node:(node.position - tyre.states.position).magnitude_squared())
            
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
            phsx.DynamicObject.add_patch(patch)
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
        phsx.DynamicObject.add_patch(rect)
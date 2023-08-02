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


class Road:
    filter_window_size_global = 0.01
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
            point_list = [np.double(line.split(",")) for line in f]
        point_list.sort(key= lambda p:p[0])
        point_list = np.array(point_list)
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
            prev_idx = next_idx = idx
            while prev_idx > 0 and (self.road.points[prev_idx] - self.position).magnitude() < self.road.filter_window_size_global:
                prev_idx -=1
            while next_idx < len(self.road.points)-1 and (self.road.points[next_idx] - self.position).magnitude() < self.road.filter_window_size_global:
                next_idx +=1
            self.tangent = (self.road.points[next_idx] -\
                                self.road.points[prev_idx]).normalized()
            self.curvature = ut.get_curvature_3point(prev = self.road.points[prev_idx],
                                                    target= self.position,
                                                    next = self.road.points[next_idx])
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
            self.start_node.section = self
            self.max_curvature_mag = np.abs(start_node.curvature)
            self.end_node:Road.Node = start_node 
            self.peak_node:Road.Node = start_node
            self.curvature_sign = np.sign(start_node.curvature)
            # find peak curvature point
            while (np.sign(self.end_node.next.curvature) == np.sign(self.end_node.curvature)) or\
                    (np.abs(self.end_node.next.curvature) < 2) :
                if np.abs(self.end_node.curvature) > self.max_curvature_mag:
                    self.peak_node = self.end_node
                    self.max_curvature_mag = np.abs(self.end_node.curvature)
                self.end_node = self.end_node.next
                self.end_node.section = self
                if self.end_node.next is None:
                    break
          
            try:
                self.fit_curvature, normal = ut.get_equivalent_circle(p1=self.start_node.position,
                                                            p0= self.peak_node.position,
                                                            p2= self.end_node.position)
                self.fit_centre = self.peak_node.position - (1/self.fit_curvature)*normal   
            except:
                self.fit_centre = None
                self.fit_curvature = None  
        def make_next(self):
            if self.next is not None:
                return self.next    
            elif self.end_node.next.next is not None:
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
                self.start_node.section = self.next
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
        def draw_fit_circle(self):
            if self.fit_centre is not None:
                fit_circle = patches.Circle(self.fit_centre,
                                            1/self.fit_curvature,
                                            fill=False, color="purple")
                phsx.DynamicObject.add_patch(fit_circle)
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
        
        initial_y = road.y[np.where(road.x > initial_x)[0][0]] + free_radius *1.1
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
        def centre_node(self):
            n = self.tyre.road.start_node
            while n.idx != self.centre_point_idx:
                n = n.next
            return n
        def draw(self):
            self.draw_terrain_circles()
            self.draw_envelop()
            pass
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
            fore_circle = patches.Wedge(center=self.fore_circle_centre,
                                        r=self.fore_circle_radius,
                                        theta1= np.rad2deg(self.centre_point_angle_f()),
                                        theta2= np.rad2deg(self.centre_point_angle_f()) + 180,
                                        fill = None,
                                        width=0.0001,
                                        lw = 2,
                                        linestyle = "--"
                                         )
            aft_circle = patches.Wedge(center = self.aft_circle_centre,
                                       r= self.aft_circle_radius,
                                       theta1 = np.rad2deg(self.centre_point_angle_f()) - 180,
                                       theta2 = np.rad2deg(self.centre_point_angle_f()),
                                       fill= None,
                                       width = 0.0001,
                                       lw =2,
                                       linestyle = "--")
            phsx.DynamicObject.add_patch(fore_circle)
            phsx.DynamicObject.add_patch(aft_circle)
        def draw_envelop(self):
            phsx.DynamicObject.plot(x= self.centre_node().position.x,
                                    y= self.centre_node().position.y,
                                    marker="x", color="green",markersize=5)
            
            x = [self.tyre.states.position.x +\
                  (self.tyre.free_radius - w)*np.cos(t) for
                    w,t in zip(self.whole_deformation_profile, self.whole_theta_profile)]
            y = [self.tyre.states.position.y +\
                  (self.tyre.free_radius - w)*np.sin(t) for
                    w,t in zip(self.whole_deformation_profile, self.whole_theta_profile)]
            phsx.DynamicObject.plot(x , y , color = "green")
        def set_equivalent_circles(self):
            fore_point = self.centre_node().section.end_node.next.position
            aft_point = self.centre_node().section.start_node.prev.position
            fore_curvature = ut.get_circle_tangent_2points(tangent=-self.normal_vector_f().cross(),
                                                           p0= self.centre_point_f(),
                                                           p1= fore_point)
            aft_curvature = ut.get_circle_tangent_2points(tangent=self.normal_vector_f().cross(),
                                                          p0 = self.centre_point_f(),
                                                          p1= aft_point)
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
        # this is a contact which takes advantage of the smart road
        # we can have multiple contacts in and the deflection from all of them is itegrated
        def __init__(self,
                     tyre, 
                     node:Road.Node,
                     ) -> None:
            self.tyre:ContinousTyre = tyre
            self.node_list: list[Road.Node] = [node]
        def get_node_angle(self , node:Road.Node):
            #no checks carried out to see if node is in the node list, BE CAREFUL
            # only call it on nodes that are in this contact
            return np.arctan2(node.position.y - self.tyre.states.position.y,
                              node.position.x - self.tyre.states.position.x)
        def get_node_normal_vector(self, node:Road.Node) -> Vector2:
            return node.position - self.tyre.states.position   
        def get_node_replacement(self, node:Road.Node):
            new_node = node
            while self.get_node_normal_vector(new_node.next).magnitude_squared() < \
                self.get_node_normal_vector(new_node).magnitude_squared():
                new_node = new_node.next
            if self.get_node_normal_vector(new_node).magnitude() < self.tyre.free_radius:
                return new_node
            return None
        def update_nodes(self):
            for i in range(len(self.node_list)):
                self.node_list[i] = self.get_node_replacement(self.node_list[i])
            self.node_list = [n for n in self.node_list if n is not None]
            return len(self.node_list) > 0
        def draw(self):
            for node in self.node_list:
                phsx.DynamicObject.plot(node.position.x, node.position.y ,
                                         marker = ".", markersize = 10, color = "purple")
                node.section.draw_fit_circle()
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
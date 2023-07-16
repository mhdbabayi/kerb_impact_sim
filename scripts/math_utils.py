import numpy as np
from euclid3 import Vector2
from scipy import io
from scipy import interpolate

def intersection(p1, p2, P1, P2):
    x0, y0 = p1
    x1, y1 = p2
    X0, Y0 = P1
    X1, Y1 = P2

    # calculate the denominator of the two equations
    den = (Y1 - Y0)*(x1 - x0) - (X1 - X0)*(y1-y0)

    # check if the lines are parallel
    if den == 0:
        return None

    # calculate the numerators of the two equations
    num1 = (X1 - X0)*(y0 - Y0) + (X0 - x0)*(Y1-Y0)
    num2 = (x1 - x0)*(y0 -Y0) + (X0 - x0)*(y1 - y0)

    # calculate the parameters the ratios of intersetion point distance from 
    # the line segment start point to line segement length
    distance_ratio_line_1 = num1 / den
    distance_ratio_line_2 = num2 / den

    # check if the intersection point is within the line segments
    if distance_ratio_line_1 >= 0 and distance_ratio_line_1 <= 1 and\
          distance_ratio_line_2 >= 0 and distance_ratio_line_2 <= 1:
        return distance_ratio_line_2

    # no intersection found
    return None
def circle_line_intersection(line_start:Vector2, line_end:Vector2, circle_centre:Vector2, radius:float):
    #v_ denotes a vector starting at line start
    v_line = line_end - line_start
    v_centre = line_start - circle_centre
    # quadratic equation constants:
    A = v_line.dot(v_line)
    B = 2*v_centre.dot(v_line)
    C = v_centre.dot(v_centre) - radius**2
    Result = [None, None]
    if (discriminant := B**2 - 4*A*C) <= 0:
        return None
    discriminant = np.sqrt(discriminant)
    t1 = (-B - discriminant)/(2*A)
    t2 = (-B + discriminant)/(2*A)
    if (t1<0) or (t1 > 1):
        if(t2<0) or (t2 > 1):
            return None
    if (t1 > 0 and t1 < 1):
        Result[0] = line_start + t1*v_line
    if(t2 > 0 and t2 < 1):
        Result[1] = line_start + t2*v_line
    return Result
def find_chord_centre(p1:Vector2, p2:Vector2, circle_centre:Vector2):
    #given two points p1 and p2, find the closest point on line segment p1-p2 to circle_centre
    # points p1 and p2 are assumed inside the circle, no checks are carried out. BE CAREFUL!
    v_line = p2 - p1
    v_centre = p1 - circle_centre
    t = -(v_line.dot(v_centre))/(v_line.magnitude_squared())
    t = min(max(t , 0),1)
    return p1 + t*v_line, t
def polar_derivative(point:Vector2, dy):
    #given point x,y and the deriavtive in cartesian coordiantes dy
    #calculates dr/d(theta) in polar coordinates
    x = point.x
    y = point.y
    r = np.sqrt(x**2 + y**2)
    dr = r*(x + y*dy)/(x*dy - y)
    return dr
def polar_second_derivative(point:Vector2, dy , ddy):

    x = point.x
    y = point.y
    r = np.sqrt(x**2 + y**2)
    dr_dtheta = polar_derivative(point, dy)
    dx_dtheta = dr_dtheta*x/r  - y
    dy_dtheta = dr_dtheta*y/r + x
    dydx_dtheta = ddy*dx_dtheta
    ddr_dtheta = ((x*dy - y)*(dx_dtheta + y*dydx_dtheta + dy*dy_dtheta) -\
                  (y*dy + x)*(dy*dx_dtheta + dydx_dtheta*x - dy_dtheta)) * r/(x*dy - y)**2\
                   +dr_dtheta*(y*dy  +x)/(x*dy - y)
    return ddr_dtheta

class BeamTyre:
    def __init__(self,
                 beta,
                 tyre_radius,
                 boundary_theta_map_file:str,
                 theta_resolution_deg = 1,
                 exp_multiplier_cutoff = 0.01) -> None:
        self.beta = beta
        self.tyre_radius = tyre_radius
        assert (exp_multiplier_cutoff < 1) and (exp_multiplier_cutoff > 0)
        max_theta = np.max((np.log(exp_multiplier_cutoff)/(-beta) , np.deg2rad(90)))
        self.theta_resolution = np.deg2rad(theta_resolution_deg)
        self.profile_theta = np.arange(start=self.theta_resolution,
                                       stop=max_theta,
                                       step=self.theta_resolution)
        self.beam_solution = lambda A, B:\
            np.exp(-self.beta*self.profile_theta)*(
            A*np.cos(beta*self.profile_theta) + B*np.sin(beta*self.profile_theta))
        self.boundary_interpolator :interpolate.RegularGridInterpolator = None
        self.terrain_radius_grid = None
        self.penetration_grid = None
        self.setup_interpolator(boundary_theta_map_file)  
    def get_profile(self, penetration, terrain_radius):
        # based on calculations in matlab file 
        theta0 = self.get_boundary_theta(penetration , terrain_radius)
        w0 = self.get_initial_displacement(penetration, terrain_radius, theta0)
        dw0 = self.get_initial_slope(penetration, terrain_radius, theta0, w0)
        return self.profile_theta + theta0, self.beam_solution(A = w0 , B = w0 + dw0/self.beta)
    def get_normal_force_integral(self,w0, dw0, theta_end = None):
        # radial force: 
        # integral exp(-beta*theta)*(A*cos(beta*theta) + B*sin(beta*theta))*cos(theta)
        # theta from 0 to theta_end. if theta_end is None, it's infinity
        beta = self.beta
        A = w0
        B = dw0/beta + w0
        integral_at_zero = beta*(2*(A + B)*beta**2 + A - B)/\
            (4*beta**4 +1)
        if theta_end is None:
            return integral_at_zero
        integral_at_theta_end = (np.exp(-beta*theta_end)/(4*beta**4+1))*\
            np.cos(beta*theta_end)*(np.sin(theta_end)*(2*B*beta**2 + A)-\
                                         beta*np.cos(theta_end)*(2*(A+B)*beta**2 + A + B))+\
            np.sin(beta*theta_end)*(beta*np.cos(theta_end)* (2*(A - B)*beta**2 + A + B) +\
                                    np.sin(theta_end)*(-2*beta**2 + A))
        return integral_at_theta_end - integral_at_zero
    def get_initial_displacement(self, penetration, terrain_radius, theta0):
        tyre_height = self.tyre_radius - penetration
        alpha = np.arcsin((1 + tyre_height/terrain_radius)*np.sin(theta0)) - theta0
        centre_distance = self.tyre_radius + terrain_radius - penetration
        return self.tyre_radius - np.sqrt(centre_distance**2 + terrain_radius**2 -\
                                        2*centre_distance*terrain_radius*np.cos(alpha))
    def get_initial_slope(self,penetration, terrain_radius, theta0, w0):
        tyre_height = self.tyre_radius - penetration
        alpha = np.arcsin((1 + tyre_height/terrain_radius)*np.sin(theta0)) - theta0
        centre_distance = self.tyre_radius + terrain_radius - penetration
        D = centre_distance
        R = terrain_radius
        dsin_alpha_dtheta = np.cos(alpha)*(D*np.cos(theta0)/(R*np.cos(alpha+theta0)) -1)
        if theta0 == 0:
            return 0
        return -R * (dsin_alpha_dtheta*np.sin(theta0) - np.cos(theta0)*np.sin(alpha))/\
            np.sin(theta0)**2
        #return -((u*D*r*np.cos(theta0)*np.sin(alpha))/(np.sqrt(1 - (u*np.sin(theta0))**2)) - 1)/\
        #    (r - w0)
    def setup_interpolator(self,file_name):
        matfile_content = io.loadmat(file_name)
        self.penetration_grid = matfile_content["penetration_num"][0]
        self.terrain_radius_grid = matfile_content["terrain_radius_num"][0]
        separation_theta = matfile_content["fitted_results"]
        self.boundary_interpolator = interpolate.RegularGridInterpolator(points= (self.penetration_grid,
                                                                                  self.terrain_radius_grid),
                                                                          values= separation_theta,
                                                                          method="linear")
    def get_boundary_theta(self, penetration, terrain_radius):
        # clip within range
        corrected_penetration = np.max((self.penetration_grid[0],
                                       np.min((self.penetration_grid[-1], penetration))))
        corrected_terrain_radius = np.max((self.terrain_radius_grid[0],
                                          np.min((self.terrain_radius_grid[-1], terrain_radius))))
        return self.boundary_interpolator((corrected_penetration,
                                           corrected_terrain_radius))
    def __call__(self, penetration , terrain_radius):
        return self.get_profile(penetration, terrain_radius)

def fit_quadratic(left_point:Vector2,
                  right_point:Vector2,
                  y0:float) -> np.polynomial.polynomial.Polynomial:
    # return quadratic numpy polynomial that passes through the point 
    # left_point, right_point and (0 , y0)
    #assert left_point.x < 0 and right_point.x > 0
    (xl, yl) = left_point
    (xr , yr)= right_point
    # TODO shouldn't have to deal with every condition separately
    if xl == xr:
        return np.polynomial.Polynomial([y0 , 0 , 0])
    if np.abs(np.rad2deg(xl)) < 1 or np.abs(np.rad2deg(xr))<1:
        return np.polynomial.polynomial.Polynomial([y0,(yr - yl)/(xr - xl)])
    p = np.polynomial.polynomial.Polynomial([
            y0,
            -(xl**2*y0 - xr**2*y0 - xl**2*yr + xr**2*yl)/(xl*xr*(xl - xr)),
            (xl*y0 - xr*y0 - xl*yr + xr*yl)/(xl*xr*(xl - xr))])
    return p
def get_quadratic_area_under(left_point:Vector2, right_point:Vector2, y0:float)-> float:
    integral_poly = fit_quadratic(left_point=left_point,
                                  right_point=right_point,
                                  y0=y0).integ()
    return integral_poly(right_point.x) - integral_poly(left_point.x)
     
def get_equivalent_circle(p1:Vector2, p0:Vector2, p2:Vector2):
    #takes three points and returns the curvature of circle that passes through the three points
    # and a vector connecting the centre of the circle to the p0
    # if the points are on a straight line (withing tolerance) return 0 for circle curvature and 
    # normal to the line for vector 
    chord_right = p2 - p0
    chord_left = p1 - p0 
    # find if surface is flat
    if (1 - np.abs(chord_right.dot(chord_left)/(chord_left.magnitude()*chord_right.magnitude()))) < 0.01:
        curvature = 0
        normal = Vector2(-chord_right.y , chord_right.x)
        return curvature , normal
    middle_point_right = p0 + chord_right/2
    middle_point_left = p0 + chord_left/2
    #lines connecting centre to middle of the chords
    radius_slope_right = Vector2(-chord_right.y , chord_right.x)
    radius_slope_left = Vector2(-chord_left.y , chord_left.x)
    # find intersection
    A = np.array([[chord_right.x , chord_right.y], 
                 [chord_left.x , chord_left.y]])
    b = np.array([middle_point_right.dot(chord_right), middle_point_left.dot(chord_left)])
    centre_x , centre_y = np.linalg.solve(A , b)
    circle_centre = Vector2(centre_x , centre_y)
    curvature = 1/((p0 - circle_centre).magnitude())
    normal = (p0 - circle_centre).normalized()
    return curvature , normal
def get_circle_tangent_2points(tangent:Vector2,
                               p0:Vector2,
                               p1:Vector2):
    normal = tangent.cross().normalized()
    curvature = np.abs(2*(p1-p0).dot(normal)/(p1-p0).magnitude_squared())
    curvature_sign = np.sign(tangent.dot(p1 - p0))
    return curvature_sign * curvature
def get_curvature_3point(prev:Vector2, target:Vector2, next:Vector2):
    # return approximate signed curvature of the circle passing all three points
    # curvature is approximated as the angle between v1 and v2 divided by arc length
    # arc approximated by 0.5(v1.mag() + v2.mag())
    # sign is positive if v2 is on the left of v1
    v1 = target - prev 
    v2 = next - target
    cross_product = v1.x*v2.y - v2.x*v1.y
    return np.arcsin(cross_product/(v1.magnitude()*v2.magnitude()))/ (0.5*(v1.magnitude() + v2.magnitude()))
def get_rotation_sign(v1:Vector2, v2:Vector2):
    return np.sign(v1.x*v2.y - v1.y*v2.x)
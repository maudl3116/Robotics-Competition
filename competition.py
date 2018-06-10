import brickpi
import time
import sys
import random
from math import *

interface = brickpi.Interface()
interface.initialize()

motors = [1, 2]
motor_sensor = [3]

interface.motorEnable(motors[0])
interface.motorEnable(motors[1])
interface.motorEnable(motor_sensor[0])

# right
motorParams1 = interface.MotorAngleControllerParameters()
motorParams1.maxRotationAcceleration = 6.0
motorParams1.maxRotationSpeed = 10.0
motorParams1.feedForwardGain = 255.0 / 25.0
motorParams1.minPWM = 25
motorParams1.pidParameters.minOutput = -255
motorParams1.pidParameters.maxOutput = 255
motorParams1.pidParameters.k_p = 500.0
motorParams1.pidParameters.k_i = 800.0
motorParams1.pidParameters.K_d = 8.0

# left
motorParams2 = interface.MotorAngleControllerParameters()
motorParams2.maxRotationAcceleration = 6.0
motorParams2.maxRotationSpeed = 10.0
motorParams2.feedForwardGain = 255.0 / 25.0
motorParams2.minPWM = 25
motorParams2.pidParameters.minOutput = -255
motorParams2.pidParameters.maxOutput = 255
motorParams2.pidParameters.k_p = 500.0
motorParams2.pidParameters.k_i = 800.0
motorParams2.pidParameters.K_d = 8.0

# up (sonar motor)
motorParams3 = interface.MotorAngleControllerParameters()
motorParams3.maxRotationAcceleration = 6.0
motorParams3.maxRotationSpeed = 3.0
motorParams3.feedForwardGain = 255.0 / 25.0
motorParams3.minPWM = 5
motorParams3.pidParameters.minOutput = -255
motorParams3.pidParameters.maxOutput = 255
motorParams3.pidParameters.k_p = 300.0
motorParams3.pidParameters.k_i = 300.0
motorParams3.pidParameters.K_d = 10.0

#touch sensors
touch_ports=[0,3]
interface.sensorEnable(touch_ports[0], brickpi.SensorType.SENSOR_TOUCH)
interface.sensorEnable(touch_ports[1], brickpi.SensorType.SENSOR_TOUCH)

# sonar setup
sonar_port = 2
interface.sensorEnable(sonar_port, brickpi.SensorType.SENSOR_ULTRASONIC);

# END PID ======================================================================

interface.setMotorAngleControllerParameters(motors[0], motorParams1)
interface.setMotorAngleControllerParameters(motors[1], motorParams2)
interface.setMotorAngleControllerParameters(motor_sensor[0], motorParams3)

turn_on_spot_constant = 0.0576320756414573
constant1 = 14.706451757812498 / 40.0

# constant10 = 0.25 * constant40

# UTILITY FUNCTIONS ####################################################



def calculate_likelihood_single_particle(x,
                                         y,
                                         theta,
                                         sonar_measurement, #  instance of SonarMeasurement
                                         the_map,
                                         the_robot):

    if sonar_measurement is None:
        return 0, False # An invalid sonar measurement gives 0 likelihood

    (m, incidence) = the_map.distance_from_particle(x, y, theta)
    
    # When incidence angle is too big, we don't update the likelihood
    # incidence_too_high = incidence > 0.7 or sonar_measurement > (100 - sonar_sensor_to_body)
    incidence_too_high = sonar_measurement.get_distance() > 150 or sonar_measurement.is_anomaly(the_robot.get_instant_position(), the_map)
    sonar_measurement_value = float(sonar_measurement.get_distance())
    m = float(m)
    
    # likelihood function 
    p_z_m = exp(-((sonar_measurement_value - m)**2.0)/(2.0 * sonar_variance)) + rubbish_k
    
    return p_z_m, incidence_too_high
    
def normalize_particles(particles):
    
    print('particles normalization')
    print('length particles ' + str(len(particles)))
    
    weights_sum = 0.0
    for (_, w) in particles:
        # print('weight for current particle is ' + str(w))
        weights_sum += w

    outparticles = []
    
    for i in range(len(particles)):
        (coords, w) = particles[i]
        w_new = w / weights_sum
        outparticles.append((coords, w_new))
        
    return outparticles

def resample_particles(particles):
    
    accumulator = 0
    
    sampling = []
    
    result = []
    
    for (coords, w) in particles:
        accumulator += w
        sampling.append((coords, accumulator))
	
    for i in range(100):
        # get random number
        randn = random.uniform(0, accumulator)
		
        for (coords, w) in sampling:
            if w >= randn:
                result.append((coords,0.01))
                break;
		
    return result

def getRandomGaussian(sigma):
    return random.gauss(0.0, sigma)

def update_particles_fwd(particles, amount):
        
   updated_particles = []

   for ((x, y, theta), w) in particles:
       e = getRandomGaussian(0.03)
       f = getRandomGaussian(0.001)
       updated_particles.append(((x + (amount + e) * cos(theta), 
                                  y + (amount + e) * sin(theta), theta + f), w))
   
   return updated_particles

# plus utilise car on utilise apply_rotation

'''def update_particles_rot_rad(particles, rads):

    updated_particles = []

    for ((x, y, theta), w) in particles:
        g = getRandomGaussian(0.005)
        updated_particles.append(((x, y, theta + rads + g), w))

    return updated_particles'''
    
# not used
'''def mean_angles(angles):
    nominator = 0
    dennominator = 0
    for a in angles:
        nominator += sin(a)
        dennominator += cos(a)
    return atan2(nominator, dennominator)'''

def is_horizontal_wall(wall):
    (ax, ay, bx, by) = wall
    return ay == by

def is_between(m, a, b):
    if a < b:
        return a <= m and m <= b
    else:
        return b <= m and m <= a

def find_minimum(xs):
    result = xs[0]
    for (x,a) in xs:
        if x < result:
            result = (x,a)
    return result

def normalize_rad_angle(r):
    return (r + (2*pi)) % (2*pi)       #useless to add 2*pi 

# ROBOT ================================================================

# A Map class containing walls
class Map:
    def __init__(self):
        self.walls = [];

    def add_wall(self,wall):
        self.walls.append(wall);

    def clear(self):
        self.walls = [];
    
    # returns a distance and also the incidence angle
    def distance_from_particle(self, x, y, theta):
        candidates = []
        for wall in self.walls:
            (ax, ay, bx, by) = wall
            m_nominator = (by-ay)*(ax-x)-(bx-ax)*(ay-y)
            m_denominator = (by-ay)*cos(theta) - (bx-ax)*sin(theta)
            try:
                m = m_nominator / m_denominator
            except:
                continue
            
            # P is the hypothetic intersection point between the particle vector and the infinite line of the wall
            px = m * cos(theta) + x
            py = m * sin(theta) + y
            
            # check if the intersection point is within the line segment of the wall
            valid = False
            if is_horizontal_wall(wall):
                valid = is_between(px, ax, bx)
            else:
                valid = is_between(py, ay, by)
            if m < 0:
                valid = False
            if valid:
                candidates.append((m, self.incidence_angle(theta, is_horizontal_wall(wall))))
        
        if len(candidates) == 0:
            print('something went terribly wrong when trying to find the distance from the wall')
            return (0,0)
        
        return find_minimum(candidates)
    
    def draw(self):
        for wall in self.walls:
            (ax, ay, bx, by) = wall
            printLine((ax, ay, 0),(bx, by, 0))
    
    # returns the expected perfect measurement from given RobotPosition
    def expected_measurement(self, robot_position, angle_offset=0):
        (dist, inc) = self.distance_from_particle(robot_position.get_x(), robot_position.get_y(), normalize_rad_angle(robot_position.get_theta() + angle_offset))
        return dist
    
    def incidence_angle(self, theta, horizontal):
        # normalize theta
        theta = normalize_rad_angle(theta)

        if horizontal:
            if theta < pi:
                return normalize_rad_angle(pi/2 - theta)
            else:
                return normalize_rad_angle(1.5*pi - theta)
        else:
            if theta < (pi/2):
                return normalize_rad_angle(theta)
            elif theta < 1.5*pi:
                return normalize_rad_angle(theta - pi)
            else:
                return normalize_rad_angle(2*pi - theta)


class RobotPosition:
    
    def __init__(self, x, y, theta):
        self.x = x
        self.y = y
        self.theta = theta
    
    def get_x(self):
        return self.x
    
    def get_y(self):
        return self.y
        
    def get_theta(self):
        return self.theta
        
    def __str__(self):
        return 'RobotPosition: {} {} {}'.format(self.x, self.y, self.theta)

class Robot:
    
    def __init__(self, init_x, init_y):
        self.n_of_particles = 100
        self.particles = [((init_x, init_y, 0.0), 1.0/self.n_of_particles) for i in range(self.n_of_particles)]   
        self.prev_motor_r, self.prev_motor_l = self.get_instant_motor_angles()
        self.step = 5 #  what is step???
        self.sonar = Sonar()
        self._map = Map()
        self._map.add_wall((0,0,0,168));        # a
        self._map.add_wall((0,168,84,168));     # b
        self._map.add_wall((84,126,84,210));    # c
        self._map.add_wall((84,210,168,210));   # d
        self._map.add_wall((168,210,168,84));   # e
        self._map.add_wall((168,84,210,84));    # f
        self._map.add_wall((210,84,210,0));     # g
        self._map.add_wall((210,0,0,0));        # h
        self.bumped = False	
        self.touch_sensor = TouchSensor()
        self.last_rotation_right = True

    def get_instant_position(self):
        x_sum = 0
        y_sum = 0
        thetas = []
        
        for ((x, y, theta), w) in self.particles:
            x_sum += x * w
            y_sum += y * w
            thetas.append((theta, w))
            
        theta_mean = self.mean_angles(thetas)
        
        return RobotPosition(x_sum, y_sum, theta_mean)
      
    def get_instant_motor_angles(self):
        motorAngles = interface.getMotorAngles(motors)
        angle_r = motorAngles[0][0]
        angle_l = motorAngles[1][0]
        return [angle_r, angle_l]
    
    # Read current state of motors, update particles
    # (similar to update_particles_fwd) 
    def apply_position(self):
        curr_angle_r, curr_angle_l = self.get_instant_motor_angles()
        
        travelled_d_r = (curr_angle_r - self.prev_motor_r) / constant1
        self.prev_motor_r = curr_angle_r

        travelled_d_l = (curr_angle_l - self.prev_motor_l) / constant1
        self.prev_motor_l = curr_angle_l

        approximated_overall_distance = -((travelled_d_r + travelled_d_l) / 2.0)
        self.particles = update_particles_fwd(self.particles, approximated_overall_distance)
        
        
    def apply_rotation(self, rad):
        updated_particles = []

        for ((x, y, theta),w) in self.particles:
            g = getRandomGaussian(0.018 * rad)
            updated_particles.append(((x, y, theta + rad + g),w))
        self.particles = updated_particles


    def rotate(self, rad):
        self.last_rotation_right = rad < 0
        
        angle = (rad * turn_on_spot_constant) * 180.0 / pi
        interface.increaseMotorAngleReferences(motors, [-angle, angle])

        while not interface.motorAngleReferencesReached(motors):
            time.sleep(0.01)

        self.apply_rotation(rad)

    # Perform Monte Carlo Localisation
    def do_MCL(self, force=False):
        # calculates likelihood of each particle
        # and update the weights
        
        # get a sonar measurement
        # :t SonarMeasurement
        sonar_measurement = self.sonar.measure_at(0) #  Measure in front
        print('Read sonar: {}'.format(sonar_measurement))
        if sonar_measurement is None:
            return
        if sonar_measurement.is_anomaly(self.get_instant_position(), self._map) and (not force):
            return
        
        likelihoods = []
        none_counter = 0
        
        for ((x,y,theta), w) in self.particles:
            likelihood, too_high = calculate_likelihood_single_particle(x, y, theta, sonar_measurement, self._map, self)
            if too_high:
                none_counter += 1
            likelihoods.append(likelihood)

        if none_counter > 50:
            print('Do not use MC')
            return
        
        outparticles = []
        for i in range(len(self.particles)):
            (coords, w) = self.particles[i]
            curr_likelihood = likelihoods[i]
            w_new = w * curr_likelihood
            outparticles.append((coords, w_new))

        outparticles = normalize_particles(outparticles)
        outparticles = resample_particles(outparticles)
        self.particles = outparticles

    
    def do_MCL_right(self):
        pass

    def normalize_rotation(self, angle):
        return ((angle + pi) % (2 * pi)) - pi

    def move_forward(self, distance_cm):
        print('Moving by {}'.format(distance_cm))
        # while moving, check for bumpers
        # if bumped, set self.bumped to True, move back 10 cm, turn 180 degrees
        # and run MCL
        constant40 = 14.706451757812498
        angle = distance_cm * constant40 / 40.0
        
        interface.increaseMotorAngleReferences(motors, [-angle, -angle])

        while not interface.motorAngleReferencesReached(motors):
            if self.touch_sensor.istouched():
                self.bumped = True
                # stop motors
                interface.setMotorRotationSpeedReferences(motors,[0,0])
                interface.setMotorPwm(motors[0],0)
                interface.setMotorPwm(motors[1],0)
                interface.increaseMotorAngleReferences(motors, [1, 1])
                while not interface.motorAngleReferencesReached(motors):
                    time.sleep(0.01)
                self.apply_position()
                break
            time.sleep(0.01)
        self.apply_position()

    def move_to_waypoint(self, input_x, input_y):
        robot_position = self.get_instant_position()
        dx = input_x - robot_position.get_x()
        dy = input_y - robot_position.get_y()
        dth = self.normalize_rotation(atan2(dy, dx) - robot_position.get_theta())
        distance = sqrt(dx**2 + dy**2)
        # rotate
        self.rotate(dth)
        # move
        self.move_forward(distance)
    
    def bump(self, expected_distance):
        driving_distance = expected_distance + 10.0
        self.move_forward(driving_distance)

    def mean_angles(self, angles):
        nominator = 0
        dennominator = 0
        for (a, w) in angles:
            nominator += sin(a)*w
            dennominator += cos(a)*w
        return atan2(nominator, dennominator)
  
    def find_obstacle(self):
        robot_position = self.get_instant_position()
        print('robot position is now {}'.format(robot_position))
        return self.sonar.find_obstacle(self._map, robot_position)
    
    def reset_sensor_head(self):
        self.sonar.rotate_to_async(0)

 
    
class SonarMeasurement:

    # Direction is relative to robot
    def __init__(self, distance, direction):
        self.distance = distance
        self.direction = direction    
    
    # returns True if this could be an anomaly (obstacle) measurement
    def is_anomaly(self, robot_position, the_map):
        if self.get_distance() > 100.0:
            return False
        expected_wall_distance = the_map.expected_measurement(robot_position, self.direction)
        if expected_wall_distance - self.get_distance() > 15.0:
            return True
        return False
    
    def get_distance(self):
        return self.distance
    
    def get_direction_relative_robot(self):
        return self.direction
    
    # robot_position is instance of RobotPosition
    def get_direction_absolute(self, robot_position):
        pass # TODO
    
    def __str__(self):
        return 'SonarMeasurement: {} {}'.format(self.distance, self.direction)
    
    
class Sonar:
    
    def __init__(self):
        self.zero_motor_angle = self.get_instant_motor_angle()
        print('Sonar zero angle is {}'.format(self.zero_motor_angle))
        self.doing_async = False
        
    def set_high_speed(self):
        # up (sonar motor)
        motorParams3 = interface.MotorAngleControllerParameters()
        motorParams3.maxRotationAcceleration = 6.0
        motorParams3.maxRotationSpeed = 3.0
        motorParams3.feedForwardGain = 255.0 / 25.0
        motorParams3.minPWM = 5
        motorParams3.pidParameters.minOutput = -255
        motorParams3.pidParameters.maxOutput = 255
        motorParams3.pidParameters.k_p = 300.0
        motorParams3.pidParameters.k_i = 300.0
        motorParams3.pidParameters.K_d = 10.0
        interface.setMotorAngleControllerParameters(motor_sensor[0], motorParams3)
        
    def set_low_speed(self):
        # up (sonar motor)
        motorParams3 = interface.MotorAngleControllerParameters()
        motorParams3.maxRotationAcceleration = 6.0
        motorParams3.maxRotationSpeed = 1.0
        motorParams3.feedForwardGain = 255.0 / 25.0
        motorParams3.minPWM = 5
        motorParams3.pidParameters.minOutput = -255
        motorParams3.pidParameters.maxOutput = 255
        motorParams3.pidParameters.k_p = 300.0
        motorParams3.pidParameters.k_i = 300.0
        motorParams3.pidParameters.K_d = 10.0
        interface.setMotorAngleControllerParameters(motor_sensor[0], motorParams3)

    # Get current direction of the sonar
    # via comparing current orientation of the motor
    # with zero_motor_angle
    # This is relative to the robot body
    # rad = 0 front of robot
    # rad = -pi/2 left
    # rad = +pi/2 right
    def get_direction(self):
        difference = self.get_instant_motor_angle() - self.zero_motor_angle
        return difference # or maybe -difference
    
    # rad is relative to robot body, as above
    def rotate_to(self, rad):
        self.rotate_to_async(rad, True)
        while not self.rotation_complete():
            time.sleep(0.01)

    def rotate_to_async(self, rad, internal=False):    #consists in reaching angle rad from the reference angle, i.e. the zeros one.
        if not internal:
            while not self.rotation_complete():
                time.sleep(0.01)
        h
        target_rotation = self.zero_motor_angle + rad
        current_angle = self.get_instant_motor_angle()
        needed_angle = target_rotation - current_angle
        print('Needed angle for sonar head is {}'.format(needed_angle))
        interface.increaseMotorAngleReferences(motor_sensor, [needed_angle,0])

    def rotation_complete(self):
        return interface.motorAngleReferencesReached(motor_sensor)
    
    # It will rotate the sonar head to specified 'rad' angle
    # and take a sonar measurement
    # rad = 0 front of robot
    # rad = -pi/2 left
    # rad = +pi/2 right
    def measure_at(self, rad):
        self.rotate_to(rad)
        try:
            reading = interface.getSensorValue(sonar_port)
            (distance, timestamp) = reading
            return SonarMeasurement(distance, self.get_direction())
        except:
            return None

    def measure_now(self):
        try:
            reading = interface.getSensorValue(sonar_port)
            (distance, timestamp) = reading
            return SonarMeasurement(distance, self.get_direction())
        except:
            return None
        
    def get_instant_motor_angle(self):
        motorAngles = interface.getMotorAngles([3,1])
        return motorAngles[0][0]
    
    def find_obstacle(self, the_map, robot_position):
        output = self.find_obstacle_helper(the_map, robot_position)
        self.set_high_speed()
        return output
    
    # the_map is instance of Map
    # robot_position is instance of RobotPosition
    def find_obstacle_helper(self, the_map, robot_position):
        # rotate sonar to -pi     # goes to the back. Does not measure, it is just reaching a place from which we can do a complete sweep
        self.rotate_to(-pi)
        self.set_low_speed()
        
        # start rotate towards pi
        self.rotate_to_async(pi)
        
        # initialize data structures
        samples = []

        print('Starting full rotation...')

        # while not reached, measure
        while not self.rotation_complete():
            # take a measurement
            sonar_measurement = self.measure_now()
            is_anomaly = sonar_measurement.is_anomaly(robot_position, the_map)
            current_angle = sonar_measurement.get_direction_relative_robot()
            samples.append((is_anomaly, current_angle))
            
            # sleep for a while
            time.sleep(0.01)
        
        print('Finished full rotation...')

        # each measurement is a tuple (is_anomaly, current_angle)
        
        # clusters
        clusters = []
        previous_cluster = []
        previous_was_cluster = True
        for s in samples:
            if s[0]:
                if not previous_was_cluster:
                    clusters.append(previous_cluster)
                    previous_cluster = []
                previous_cluster.append(s)
            previous_was_cluster = s[0]
        
        i = 0
        
        # find the largest cluster
        largest_index = -1
        largest_size = None
        i = 0
        for i in range(len(clusters)):
            c = clusters[i]
            if largest_size is None or len(c) > largest_size:
                if len(c) > 0 and (c[len(c)-1][1] - c[0][1] > (pi/2)):
                    continue
                largest_index = i
                largest_size = len(c)
        try:
            largest_cluster = clusters[largest_index]
        except:
            return None
        
        # guess which is an obstacle
        print('largest cluster has {}'.format(largest_cluster))
        
        angles = []
        for s in largest_cluster:
            if s[0]:
                angles.append(s[1])
        
        if len(angles) == 0:
            return None
        
        # calculate mean angle of cluster
        mean_sum = 0.0
        for a in angles:
            mean_sum += a
        return mean_sum / len(angles)
    

# TOUCH SENSOR #########################################################

class TouchSensor():
    
    def __init__(self):
        pass
    
    def touched(self):
        resultL=interface.getSensorValue(touch_ports[0])
        resultR=interface.getSensorValue(touch_ports[1])
        touchedL=resultL[0]
        touchedR=resultR[0]
        return [touchedL,touchedR]

    def istouched(self):
        touch_result = self.touched()
        return touch_result[0] or touch_result[1]	



# MONTE CARLO PARAMETERS =======================================================

sonar_variance = 2.0 # Variance of sonar measurements
rubbish_k = 0.0005 # Up-lift of Gaussian to allow small probability of completely rubbish measurements from sonar
sonar_sensor_to_body = 2.0 # distance from sonar to robot body centre offset


# Printing utilities ###################################################

#def scale_2d_point(x, y):
 #   return (x*4+100, 500-(y*4+100))
def scale_2d_point(x, y):
    return (x*3+100, 700-(y*3+50))

    
def printParticles(particles):
    # ((x, y, theta), w)
    newarray = []
    for ((x, y, theta), w) in particles:
        (nx, ny) = scale_2d_point(x, y)
        newparticle = (nx, ny, theta)
        newarray.append(newparticle)
    print 'drawParticles:' + str(newarray)
    time.sleep(0.1)
    
def printLine(current, new):
    (cx, cy, _) = current
    (nx, ny, _) = new
    (xc, yc) = scale_2d_point(cx, cy)
    (xn, yn) = scale_2d_point(nx, ny)
    print 'drawLine:' + str((xc, yc, xn, yn))


# MAIN ALGORITHM #######################################################

# TODO put sensible stuff here


def analyze_area(the_robot, area_coords, exit_coords):
    print('Starting to analyze area')
    
    sensor_angle_offset = 0.9
    
    the_robot.bumped = False
    
    centre_x, centre_y = area_coords
    the_robot.move_to_waypoint(centre_x, centre_y)
    the_robot.do_MCL()
    
    if the_robot.bumped:
        the_robot.move_forward(-10.0)
        if the_robot.last_rotation_right:
            the_robot.rotate(pi/2)
        else:
            the_robot.rotate(-pi/2)
        the_robot.do_MCL()
        return # if bumps early, return directly

    obstacle_direction = the_robot.find_obstacle()    # fait le sweep, et renvoie la direction dans laquelle a vu la bouteille
    if obstacle_direction is None:
        exit_x, exit_y = exit_coords
        the_robot.move_to_waypoint(exit_x, exit_y)
        the_robot.do_MCL()
        print('Area skipped, could not find obstacle')
        return
    print('Obstacle found at {}'.format(obstacle_direction))
    
    # start rotate sensor head to front position
    the_robot.reset_sensor_head()
    
    # rotate by obstacle direction
    the_robot.rotate(obstacle_direction * sensor_angle_offset)   # la direction a deja bien ete calculee, pour qu'on aille vers la bouteille
    the_robot.bump(70)
    print('Bumped into obstacle!')
    the_robot.move_forward(-10.0)
    
    exit_x, exit_y = exit_coords
    the_robot.move_to_waypoint(exit_x, exit_y)
    the_robot.do_MCL()
    
    print('Area complete')

def main_competition():
    area_A_coords = (160.0, 40.0)
    area_B_coords = (126.0, 147.0)
    area_C_coords = (42.0, 126.0)
    
    exit_A_coords = (126.0, 40.0)
    exit_B_coords = (126.0, 84.0)
    exit_C_coords = (42.0, 84.0)
    
    the_robot = Robot(84.0, 30.0)
    analyze_area(the_robot, area_A_coords, exit_A_coords)
    analyze_area(the_robot, area_B_coords, exit_B_coords)
    analyze_area(the_robot, area_C_coords, exit_C_coords)
    
    the_robot.do_MCL(True)
    the_robot.move_to_waypoint(42.0, 42.0)
    the_robot.do_MCL(True)
    the_robot.rotate(-pi/2)
    the_robot.do_MCL(True)
    the_robot.move_to_waypoint(84.0, 30.0)

def test_face_wall():
    the_robot = Robot(0.0, 0.0)
    the_sonar = Sonar()
    theta = (the_robot.get_instant_position()).get_theta()
    theta_90 = the_robot.normalize_rotation(-pi/2 - theta)
    the_robot.do_MCL_head(False,theta_90)
    theta_180 = the_robot.normalize_rotation(-pi - theta_90)
    the_robot.do_MCL_head(False,-pi)#theta_180)
    #the_sonar.rotate_to(theta_180)
    
def main_forward_bump_test():
    the_robot = Robot(0.0, 0.0)
    the_robot.move_forward(40.0)
    
def main_sonar():
    the_sonar = Sonar()
    the_sonar.rotate_to(-pi)
    the_sonar.rotate_to(pi)

def main_find_obstacle():
    the_robot = Robot(84.0, 30.0)
    the_robot.move_to_waypoint(160.0, 40.0)
    obstacle_direction = the_robot.find_obstacle()
    print('Obstacle found at {}'.format(obstacle_direction))
    # rotate by obstacle direction
    the_robot.rotate(-obstacle_direction)
    the_robot.move_forward(40.0)
    

def main_rotation():
    the_robot = Robot(0.0, 0.0)
    the_robot.rotate(pi)

def main_bump():
    the_robot = Robot(84.0, 30.0)
    the_robot.move_forward(70)
    the_robot.move_forward(-20)
  
def main_touch():
    the_robot = Robot(84.0, 30.0)
    while True:
        print(the_robot.touch_sensor.istouched())
        time.sleep(0.1)
    
# main_forward_bump_test()
main_competition()
#test_face_wall()
# main_sonar()
# main_find_obstacle()
# main_rotation()
# main_bump()
# main_touch()

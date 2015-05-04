import time, sys, math, csv
from math import log10, sin, cos, floor, asin, pi, sqrt, fabs, exp

# 2D Buoyancy script
# A rod with length L has a odd number set of test points, the middle point representing the origin. 
# we will assume the rod is made of pine and the liquid is water therefore having a density ratio of 2
class Buoyancy:
    def __init__(self, numpoints, waterHeight, initialVelocity, rodMass, rodLength, rotation):
        self._debugging = False
        if self._debugging:
            self._timeInterval = 0.1
            self._timeEnd = 10.0
        else:
            self._timeInterval = .001
            self._timeEnd = 30.0

        self._fluidDensity = 1.0
        self._objectDensity = 0.5

        # For Drag force, b value -> (density of medium)*(cross sectional area)
        self._waterLinearDampingFactor = 1.0
        self._waterAngularDampingFactor = 10.0
        self._airLinearDampingFactor = 0.1
        self._airAngularDampingFactor = 1.0

        self._rounding = True
        self._sigfig = 4
        self._waterheight = -500.0
        self._floor = -200.0
        self._mInertia = (rodMass*pow(rodLength,2.0))/12.0
        angle_rad = self.convertDegToRad(rotation)
        points = self.generate_points(numpoints, rodLength)
        points = self.TransformPoints(points, angle_rad)
        data = self.simulate_buoyancy(points, rodMass, rodLength, initialVelocity, rotation)
        if self._rounding:
            data = self.roundData(data)

        if self._debugging:
            self.evaluate(data)
        else:
            self.writeCSV(data)
            self.evaluate(data)

    def convertRadToDeg(self, radians):
        degrees = radians * 57.2957795
        return degrees

    def convertDegToRad(self, degrees):
        radians = degrees/57.2957795
        return radians

    def generate_points(self, numpoints, rodLength):
        points = []
        if numpoints < 3:
            print 'ERROR: invalid number of points.'
            sys.exit(0),'\n'

        if numpoints % 2 == 0:
            #even number, we want odd, add +1, notify user
            numpoints = numpoints + 1
            print 'WARNING: Even number of points was given, an additional point has been added.'

        i = 0.0 
        #generate points along the x-axis
        while i < numpoints:
            x = -rodLength/2.0 + ((rodLength*i)/(float(numpoints)-1.0))
            #create a point that lies on the x-axis
            point = [x,0]
            points.append(point)
            i+=1.0
            
        return points

    def simulate_buoyancy(self, points, rodMass, rodLength, initialVelocity, angle):
        # Initial parameters/instantiating data lists
        t = 0.0
        angular_velocity = 0.0
        velocity = initialVelocity
        numPoints = len(points)
        density_ratio = self._fluidDensity/self._objectDensity
        # Initial acceleration is downward gravity, if above water
        if points[numPoints/2][1] > self._waterheight:
            acceleration = -9.8
        else:
            acceleration = density_ratio*9.8
        print acceleration
        data = {}
        times = []
        newpoints_sum = []
        torques = [] 
        forces = []
        angular_accelerations = []
        angular_velocitys = [] 
        linear_accelerations = []
        linear_velocitys = []
        angles = []
        origin_pointsX = []
        origin_pointsY = []
        # debugging purposes
        segments = []
        segment_sums = []
        isStraight_list = []
        
        # convert angle from x axis into rads, math.cos/sin only take rad in puts
        # theta is the initial angle of rotation
        theta = self.convertDegToRad(angle)
        
        while t < self._timeEnd:

            # reset torque and force sum for every time iteration
            torque_sum = 0.0
            force_sum = 0.0

            #Debugging Purposes
            segment_sum = 0.0
            segment = []
            isStraight = []
            i = 0

            newpoints = []
            for point in points:
                # Drag force is initially calculated in the +y direction, assumign movement is in the -y direction
                Fg = ((rodMass/numPoints)*9.8)
                Fb = ((rodMass/numPoints)*density_ratio*9.8)
                Flinear_drag = pow(velocity,2.0)/float(numPoints)
                if velocity > 0.0:
                    Flinear_drag = -Flinear_drag
                y = velocity*self._timeInterval + pow(self._timeInterval,2.0) *(0.5*acceleration) + point[1]
                
                if y < self._waterheight:
                    # Object's Net Force = Buoyancy Force + Drag Force - Force of Gravity
                    Flinear_drag = Flinear_drag * self._waterLinearDampingFactor
                    Fnet = Fb + Flinear_drag - Fg
                    angDrag_coe = self._waterAngularDampingFactor
                else:
                    # Net Force = - Force of Gravity on object
                    Flinear_drag = Flinear_drag * self._airLinearDampingFactor
                    Fnet = -Fg + Flinear_drag
                    angDrag_coe = self._airAngularDampingFactor
                newpoint = [point[0],y]
                newpoints.append(newpoint)
                #get center point of rod
                rod_origin = points[numPoints/2]

                # Calculate displacement vector, r
                r_x = point[0] - rod_origin[0]
                r_y = point[1] - rod_origin[1]
                r = sqrt(pow(r_x,2.0) + pow(r_y,2.0))

                # Calculate Tangential Velocity and angular drag force, w*r = V, Fd = b2*V^2 
                V_tan = angular_velocity * r
                Fangular_drag = angDrag_coe*pow(V_tan,2.0)

                # Torque = r x F
                torque = r_x*Fnet
                torque_drag = r * Fangular_drag
                if angular_velocity > 0.0:
                    torque_drag = -torque_drag
                # Drag Torque points opposite of angular velocity
                torque_net = torque + torque_drag
                torque_sum += torque_net
                force_sum += Fnet

                #For Debugging, get length of segment
                if i < numPoints - 1:
                    next_point = points[i+1]
                    x_len = fabs(point[0] - next_point[0])
                    y_len = fabs(point[1] - next_point[1])
                    length = sqrt(pow(x_len,2.0) + pow(y_len,2.0))
                    segment.append(length)
                    segment_sum += length
                i+=1

                # end of for loop


            #Update linear and angular acceleration/velocity and the angle of rotation of the rod
            acceleration = force_sum/rodMass
            velocity = self.UpdateLinearVelocity(velocity, self._timeInterval, acceleration)
            angular_acceleration = torque_sum/self._mInertia
            angular_velocity = self.UpdateAngularVelocity(angular_velocity, self._timeInterval, angular_acceleration)
            #We might want to calculate theta before angular velocity, caution
            theta = self.calculateAngle(angular_velocity, self._timeInterval, angular_acceleration) 
            newpoints = self.TransformPoints(newpoints, theta)
            isStraight = self.isStraight(points, rodLength)
            points = newpoints
            angle_change = self.convertRadToDeg(theta)
            angle = angle + angle_change
                
            times.append(t)
            segments.append(segment)
            segment_sums.append(segment_sum)
            isStraight_list.append(isStraight)
            newpoints_sum.append(newpoints)
            torques.append(torque_sum)
            forces.append(force_sum)
            angular_accelerations.append(angular_acceleration)
            angular_velocitys.append(angular_velocity)
            linear_accelerations.append(acceleration)
            linear_velocitys.append(velocity)
            angles.append(angle)
            origin_pointsX.append(points[numPoints/2][0])
            origin_pointsY.append(points[numPoints/2][1])

            data = {
            'time':times,
            'point locations':newpoints_sum,
            'segments':segments,
            'is straight/length difference':isStraight_list,
            'rod length':segment_sums,
            'torque': torques,
            'force':forces,
            'angular velocity':angular_velocitys,
            'angular acceleration':angular_accelerations,
            'linear acceleration':linear_accelerations,
            'linear velocity':linear_velocitys, 
            'angle':angles,
            'originX':origin_pointsX,
            'originY':origin_pointsY,
            }

            t += self._timeInterval

        return data

    def calculateAngle(self, angular_velocity, time, angular_acceleration):
        return ((angular_velocity*time) + (0.5*angular_acceleration*pow(time,2.0))) 

    def UpdateAngularVelocity(self, angular_velocity, time, angular_acceleration):
        return (angular_velocity + angular_acceleration*time)

    def UpdateLinearVelocity(self, velocity, time, acceleration):
        return (velocity + acceleration*time)

    def isStraight(self, points, rodLength):
        # This checks if the points line up in a straight line by taking the end points and calculating the length and referencing it to the given rod length
        data = []
        numpoints = len(points)
        first_point = points[0]
        last_point = points[numpoints -1]
        x_len = fabs(first_point[0] - last_point[0])
        y_len = fabs(first_point[1] - last_point[1])
        length = sqrt(pow(x_len,2.0) + pow(y_len,2.0))
        difference = fabs(rodLength - length)
        if self._rounding:
            length = round(length,self._sigfig)
            difference = round(difference,self._sigfig)
        if length == rodLength:
            isStraight = True
        else:
            isStraight = False

        data.append(isStraight)
        data.append(difference)
        return data

    def roundData(self, data):
        times = data['time']
        angular_velocitys = data['angular velocity']
        angular_accelerations = data['angular acceleration']
        linear_accelerations = data['linear acceleration']
        linear_velocitys = data['linear velocity']
        locations = data['point locations']
        torques = data['torque']
        forces = data['force']
        angles = data['angle']
        segments = data['segments']
        segment_sums = data['rod length']
        isStraight_list = data['is straight/length difference']
        origin_pointsX = data['originX']
        origin_pointsY = data['originY']
        i = 0
        while i < len(times):

            times[i] = round(times[i],self._sigfig)
            angular_velocitys[i] = round(angular_velocitys[i],self._sigfig)
            angular_accelerations[i] = round(angular_accelerations[i],self._sigfig)
            linear_accelerations[i] = round(linear_accelerations[i],self._sigfig)
            linear_velocitys[i] = round(linear_velocitys[i],self._sigfig)
            torques[i] = round(torques[i],self._sigfig)
            forces[i] = round(forces[i],self._sigfig)
            angles[i] = round(angles[i],self._sigfig)
            segment_sums[i] = round(segment_sums[i],self._sigfig)
            points = locations[i]
            seg = segments[i]
            isStraight = isStraight_list [i]

            isStraight[1] = round(isStraight[1],self._sigfig)

            if isStraight[1] == 0.0:
                isStraight[1] = 'No difference'
            else:
                isStraight[1] = repr(difference) + ' difference'

            # 5 loops
            for point in points:
                point[0] = round(point[0],self._sigfig)
                point[1] = round(point[1],self._sigfig)
            # 4 loops
            for s in seg:
                s = round(s,self._sigfig)

            locations[i] = points
            segments[i] = seg
            isStraight_list[i] = isStraight

            i+=1

        new_data = {'time':times,
        'point locations':locations,
        'segments':segments,
        'is straight/length difference':isStraight_list,
        'rod length':segment_sums,
        'torque': torques,
        'force':forces,
        'angular velocity':angular_velocitys,
        'angular acceleration':angular_accelerations,
        'linear acceleration':linear_accelerations,
        'linear velocity':linear_velocitys, 
        'angle':angles,
        'originX':origin_pointsX,
        'originY':origin_pointsY}
        return new_data


    def TransformPoints(self, points, theta):
        # Preforms a 2D counter clockwise rotation with respect to the angle from the x axis
        
        #Get the center point of the rod
        rod_origin = points[len(points)/2]
        #Get the distance from the center point of the rod to the origin for each direction
        x_dist = rod_origin[0]
        y_dist = rod_origin[1]
        
        # Rotational Matrix becomes:
        # x' = xcos - ysin
        # y' = xsin + ycos

        for point in points:
            #Translate point to the origin
            x = point[0]
            y = point[1]
            x = x - x_dist
            y = y - y_dist

            #Make rotation on point with respect to origin using the rotation matrix
            x_prime = x*cos(theta) - y*sin(theta)
            y_prime = x*sin(theta) + y*cos(theta)

            #Return point to its original origin
            x_transformed = x_prime + x_dist
            y_transformed = y_prime + y_dist

            point[0] = x_transformed
            point[1] = y_transformed

        return points

    def writeCSV(self, data):
        with open('2DBuoyancy.csv', 'wb') as csvfile:
            writer = csv.DictWriter(csvfile, data.keys())
            writer.writeheader()
            times = data['time']
            angular_velocitys = data['angular velocity']
            angular_accelerations = data['angular acceleration']
            linear_accelerations = data['linear acceleration']
            linear_velocitys = data['linear velocity']
            locations = data['point locations']
            torques = data['torque']
            forces = data['force']
            angles = data['angle']
            segments = data['segments']
            segment_sums = data['rod length']
            isStraight_list = data['is straight/length difference']
            origin_pointsX = data['originX']
            origin_pointsY = data['originY']

            i = 0
            while i < len(times):
                row = {'force':forces[i],
                        'torque':torques[i],
                        'point locations':locations[i],
                        'segments':segments[i],
                        'rod length':segment_sums[i],
                        'is straight/length difference':isStraight_list[i],
                        'angular velocity':angular_velocitys[i],
                        'angular acceleration':angular_accelerations[i],
                        'linear velocity':linear_velocitys[i],
                        'linear acceleration':linear_accelerations[i],
                        'time':times[i],
                        'angle':angles[i],
                        'originX':origin_pointsX[i],
                        'originY':origin_pointsY[i]
                        }
                writer.writerow(row)
                i+=1

    def evaluate(self, data):
        times = data['time']
        angular_velocitys = data['angular velocity']
        angular_accelerations = data['angular acceleration']
        linear_accelerations = data['linear acceleration']
        linear_velocitys = data['linear velocity']
        locations = data['point locations']
        torques = data['torque']
        forces = data['force']
        angles = data['angle']
        segments = data['segments']
        segment_sums = data['rod length']
        isStraight_list = data['is straight/length difference']

        # Find the terminal velocity when the rod hits the water, compare it to the max velocity. 
        # The idea is to use the damping factor to make the terminal velocity be the max velocity
        t_water = self.getTimeHitsWater(locations)
        max_velocity = self.getMaxVelocity(linear_velocitys)
        terminal_velocity = linear_velocitys[int(t_water/self._timeInterval) - 1]

        prompt = 'Terminal Velocity: ' + repr(terminal_velocity) + ', Max Velocity: ' + repr(max_velocity)
        if fabs(terminal_velocity) < fabs(max_velocity):
            prompt = prompt + ' Consider inreasing damping.' 
        print prompt
    
    def getTimeHitsWater(self, locations):
        t = 0
        while t < len(locations):
            points = locations[t]
            for point in points:
                if point[1] < self._waterheight:
                    return t*self._timeInterval
            t+=1

    def getMaxVelocity(self,linear_velocitys):
        max_velocity = 0.0
        for velocity in linear_velocitys:
            if fabs(velocity) > fabs(max_velocity):
                max_velocity = velocity
        return max_velocity


static_test = Buoyancy(5, -50.0, 0.0, 10.0, 10.0, -36.87)

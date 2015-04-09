import time, sys, math, csv
from math import log10, sin, cos, floor, asin, pi, sqrt, fabs

# 2D Buoyancy script
# A rod with length L has a odd number set of test points, the middle point representing the origin. 
# we will assume the rod is made of pine and the liquid is water therefore having a density ratio of 2
class Buoyancy:
    def __init__(self, numpoints, waterHeight, initialVelocity, rodMass, rodLength, rotation):
        self._debugging = True
        if self._debugging:
            self._timeInterval = 0.1
        else:
            self._timeInterval = .0005
        self._timeEnd = 4.0
        self._fluidDensity = 1.0
        self._objectDensity = 0.5
        self._rounding = True
        self._sigfig = 3
        self._waterheight = -50.0
        angle_rad = self.convertDegToRad(rotation)
        points = self.generate_points(numpoints, rodLength)
        points = self.TransformPoints(points, angle_rad)
        data = self.simulate_buoyancy(points, rodMass, rodLength, initialVelocity, rotation)
        if self._debugging:
            self.promptData(data)
        else:
            self.writeCSV(data)

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
        #instantiate a float numpoints
        fnumpoints = numpoints * 1.0
        #generate points along the x-axis
        while i < numpoints:
            x = -rodLength/2.0 + ((rodLength*i)/(fnumpoints-1.0))
            #create a point that lies on the x-axis
            if self._rounding:
                x = round(x,self._sigfig)
            point = [x,0]
            points.append(point)
            i+=1.0
            
        return points

    def simulate_buoyancy(self, points, rodMass, rodLength, initialVelocity, angle):
        # Initial parameters/instantiating data lists
        t = 0.0
        ti = 0.0
        angular_velocity = 0.0
        velocity = initialVelocity
        free_fall = True
        numPoints = len(points)
        yi = self.getInitHeight(points)
        density_ratio = self._fluidDensity/self._objectDensity
        # Initial acceleration is downward gravity
        acceleration = -9.8
        data = {}
        times = []
        newpoints_sum = []
        torques = [] 
        forces = []
        angular_accelerations = [] 
        linear_accelerations = []
        angles = []
        
        # convert angle from x axis into rads, math.cos/sin only take rad in puts
        # theta is the initial angle of rotation
        theta = self.convertDegToRad(angle)
        
        while t < self._timeEnd:

            # change in time from freefall
            deltaT = t - ti

            # reset torque and force sum for every time iteration
            torque_sum = 0.0
            force_sum = 0.0
            
            # j keeps track on which point we are on in the for the point array loop
            j = 0

            newpoints = []

            for point in points:
                #height transform, Force is in the y-direction
                #WARNING: This code does not update velocity
                y = velocity*t + pow(t,2.0) *(0.5*acceleration) + yi[j] 
                
                if y < self._waterheight:
                    # Net Force = Buoyancy Force on object - Force of Gravity on object
                    Fnet = ((rodMass/numPoints)*density_ratio*9.8) - ((rodMass/numPoints)*9.8)     
                else:
                    # Net Force = - Force of Gravity on object
                    Fnet = -(rodMass/numPoints)*9.8

                if self._rounding:
                    y = round(y,self._sigfig)
                newpoint = [point[0],y]
                newpoints.append(newpoint)
                #get center point of rod
                rod_origin = points[numPoints/2]
                r = point[0] - rod_origin[0]
                torque = r*Fnet
                torque_sum += torque
                force_sum += Fnet
                j += 1

                # end of for loop

            #Update acceleration of the rod
            acceleration_temp = force_sum/rodMass
            if acceleration_temp != acceleration:
                #acceleration has changed 
                acceleration = acceleration_temp
                ti = t
                yi = self.getInitHeight(points)
            # Calculate the angular acceleration using torque divided by the moment of inertia for that object
            angular_acceleration = torque_sum/((rodMass*pow(rodLength,2.0))/12.0)
            # Update theta using the last angular_velocity, for time we want to use the change in time which is the time iteration variable
            theta = self.calculateAngle(angular_velocity, self._timeInterval, angular_acceleration) 
            # Update the angular_velocity if angular acceleration is changing, for time we want to use the change in time
            angular_velocity = self.UpdateAngularVelocity(angular_velocity, self._timeInterval, angular_acceleration)
            newpoints = self.TransformPoints(newpoints, theta)
            Valid = self.CheckLength(points, rodLength)
            print Valid
            points = newpoints
            angle_change = self.convertRadToDeg(theta)
            angle = angle - angle_change

            if self._rounding:
                torque_sum = round(torque_sum,self._sigfig)
                force_sum = round(force_sum,self._sigfig)
                angle = round(angle,self._sigfig)
                angular_acceleration = round(angular_acceleration,self._sigfig) 
                acceleration = round(acceleration,self._sigfig)
                angle = round(angle,self._sigfig)
                
            times.append(t)
            newpoints_sum.append(newpoints)
            torques.append(torque_sum)
            forces.append(force_sum)
            angular_accelerations.append(angular_acceleration)
            linear_accelerations.append(acceleration)
            angles.append(angle)

            data = {'time':times,'point locations':newpoints_sum,'torque': torques,'force':forces,'angular acceleration':angular_accelerations,'linear acceleration':linear_accelerations, 'angle':angles}

            t += self._timeInterval

        return data

    def getInitHeight(self, points):
        # Initial height of points.
        yi = []
        for point in points:
            yi.append(point[1])
        return yi

    def calculateAngle(self, angular_velocity, time, angular_acceleration):
        return ((angular_velocity*time) + (0.5*angular_acceleration*pow(time,2.0))) 

    def UpdateAngularVelocity(self, angular_velocity, time, angular_acceleration):
        return (angular_velocity + angular_acceleration*time)

    def CheckLength(self, points, rodLength):
        numpoints = len(points)
        length = sqrt(pow(fabs(points[0][0]) + fabs(points[numpoints-1][0]),2.0) + pow(fabs(points[0][1]) + fabs(points[numpoints-1][1]),2.0))
        if length == rodLength:
            return True
        else:
            return False
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

            if self._rounding:
                x_transformed = round(x_transformed,self._sigfig)
                y_transformed = round(y_transformed,self._sigfig)
            point[0] = x_transformed
            point[1] = y_transformed

        return points

    def writeCSV(self, data):
        with open('2DBuoyancy.csv', 'wb') as csvfile:
            writer = csv.DictWriter(csvfile, data.keys())
            writer.writeheader()
            times = data['time']
            angular_accelerations = data['angular acceleration']
            linear_accelerations = data['linear acceleration']
            locations = data['point locations']
            torques = data['torque']
            forces = data['force']
            angles = data['angle']
            i = 0
            while i < len(times):
                row = {'force':forces[i],
                        'torque':torques[i],
                        'point locations':locations[i],
                        'angular acceleration':angular_accelerations[i],
                        'linear acceleration':linear_accelerations[i],
                        'time':times[i],
                        'angle':angles[i],
                        }
                writer.writerow(row)
                i+=1

    def promptData(self, data):
        times = data['time']
        angular_accelerations = data['angular acceleration']
        linear_accelerations = data['linear acceleration']
        locations = data['point locations']
        torques = data['torque']
        forces = data['force']
        angles = data['angle']
        i = 0
        while i < len(times):
            #print ('time:',times[i],'points:',locations[i],'angular acceleration:',angular_accelerations[i],'linear acceleration:',linear_accelerations[i],'torque:',torques[i],'force:',forces[i],'angle:',angles[i])
            i+=1

static_test = Buoyancy(5, -50.0, 0.0, 10.0, 10.0, -36.87)

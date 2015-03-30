import time, sys, math, csv
from math import log10, sin, cos, floor

# 2D Buoyancy script
# A rod with length L has a odd number set of test points, the middle point representing the origin. 
# we will assume the rod is made of pine and the liquid is water therefore having a density ratio of 2
class Buoyancy:
    def __init__(self, numpoints, waterHeight, initialVelocity, rodMass, rodLength, rotation):
        angle_rad = self.convertDegToRad(rotation)
        self._rounding = True
        self._sigfig = 3
        points = self.generate_points(numpoints, rodLength, angle_rad)
        data = self.simulate_buoyancy(points, rodMass, rodLength, initialVelocity, waterHeight, rotation)
        self.writeCSV(data)

    def convertRadToDeg(self, radians):
        degrees = radians * 57.2957795
        return degrees

    def convertDegToRad(self, degrees):
        radians = degrees/58.2957795
        return radians

    def generate_points(self, numpoints, rodLength, rotation):
        points = []
        if numpoints < 3:
            print 'ERROR: invalid number of points.'
            sys.exit(0),'\n'

        if numpoints % 2 == 0:
            #even number, we want odd, add +1, notify user
            numpoints = numpoints + 1
            print 'WARNING: Even number of points was given, an additional point has been added.'

        x_rot = cos(rotation)
        y_rot = sin(rotation)

        i = 0.0 
        while i < numpoints:
            nump = numpoints * 1.0
            increment = (rodLength*i)/((numpoints*1.0)-1.0)
            x_coe = -rodLength/2.0 + increment
            y_coe = rodLength/2.0 - increment
            if self._rounding:
                x = round(x_coe * x_rot, self._sigfig)
                y = round(y_coe * y_rot, self._sigfig)
            else:
                x = x_coe * x_rot
                y = y_coe * y_rot
            point = [x,y]
            points.append(point)
            i+=1.0

        return points

    def simulate_buoyancy(self, points, rodMass, rodLength, initialVelocity, waterHeight, rotation):
        t = 0.0
        numPoints = len(points)
        yi = self.getInitHeight(points)
        data = {}
        times = []
        newpoints_sum = []
        torques = [] 
        forces = []
        angularmomentums = [] 
        while t < 4.0:
            torque_sum = 0.0
            force_sum = 0.0
            # j keeps track on which point we are on
            j = 0

            newpoints = []

            for point in points:
                #height transform, Force is in the y-direction
                if self._rounding:
                    y = round(initialVelocity*t + pow(t,2.0) *(0.5*-9.8) + yi[j],self._sigfig) 
                else:
                    y = initialVelocity*t + pow(t,2.0) *(0.5*-9.8) + yi[j] 
                if y < waterHeight:
                    #the density ratio is 1/0.5 = 2 
                    if self._rounding:
                        Fnet = round(((rodMass/numPoints)*2.0*9.8) - ((rodMass/numPoints)*9.8),self._sigfig) 
                    else:
                        Fnet = ((rodMass/numPoints)*2.0*9.8) - ((rodMass/numPoints)*9.8)     
                    #point[1] = -50.0
                    #newpoint = [point[0],-50.0]
                    point[1] = y
                    newpoint = [point[0],y]
                else:
                    if self._rounding:
                        Fnet = round(-(rodMass/numPoints)*9.8,self._sigfig)
                    else:
                        Fnet = -(rodMass/numPoints)*9.8
                    point[1] = y
                    newpoint = [point[0],y]
                newpoints.append(newpoint)
                #get center point of rod
                rod_origin = points[numPoints/2]
                # Torque = r * F sin(theta)
                # Torque in the y direction
                # convert angle from x axis into rads, math.cos/sin only take rad in puts
                theta = self.convertDegToRad(rotation)
                r = point[1] - rod_origin[1]
                if self._rounding:
                    torque = round(r*Fnet*sin(theta),self._sigfig)
                else:
                    torque = r*Fnet*sin(theta)
                torque_sum += torque
                force_sum += Fnet
                j += 1
            if torque_sum < 0.1:
                torque_sum = 0.0
            if self._rounding:
                angularmomentum = round(torque_sum/((rodMass*pow(rodLength,2.0))/12.0),self._sigfig) 
            else:
                angularmomentum = torque_sum/((rodMass*pow(rodLength,2.0))/12.0)
            oldpoints = points
            times.append(t)
            newpoints_sum.append(newpoints)
            # torque is being rounded to wierd values this clamps it, figure out the problem.
            torques.append(torque_sum)
            forces.append(force_sum)
            angularmomentums.append(angularmomentum)

            data = {'time':times,'point locations':newpoints_sum,'torque': torques,'force':forces,'angular momentum':angularmomentums, 'angle':rotation}

            t += 0.0005

        return data

    def getInitHeight(self, points):
        # Initial height of points.
        yi = []
        for point in points:
            yi.append(point[1])
        return yi

    def writeCSV(self, data):
        with open('2DBuoyancy.csv', 'wb') as csvfile:
            writer = csv.DictWriter(csvfile, data.keys())
            writer.writeheader()
            times = data['time']
            angularmomentums = data['angular momentum']
            locations = data['point locations']
            torques = data['torque']
            forces = data['force']
            angle = data['angle']
            i = 0
            while i < len(times):
                row = {'force':forces[i],
                        'torque':torques[i],
                        'point locations':locations[i],
                        'angular momentum':angularmomentums[i],
                        'time':times[i],
                        'angle':angle,
                        }
                writer.writerow(row)
                i+=1

static_test = Buoyancy(5, -50.0, 0.0, 10.0, 10.0, 143.13)

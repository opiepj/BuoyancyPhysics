""" 
    3DBuoyancy.py
    Patrick Opie 
    Simple script that simulates buoyancy of a specified rod in 3D space. The data is outputted to a csv file.
"""
import csv, vector, math
from vector import Vector
from math import sin, cos, pi, sqrt, fabs, exp

class Buoyancy:
    def __init__(self, numpoints, water_plane, init_velocity, obj_mass, obj_length, init_rot):
        self._time_interval = 0.001
        self._time_end = 30.0
        self._fluid_density = 1.0
        self._obj_density = 0.5

        self._rodmass = obj_mass
        self._rodlength = obj_length

        self._water_linear_damping = 1.0
        self._water_angular_damping = 10.0
        self._air_linear_damping  = 0.1
        self._air_angular_damping = 1.0

        self._water_plane_z = -200 #Height of water plane (z-dir)

        self._moment_inertia = (obj_mass*pow(length,2.0))/12.0
        points = self.generatepoints(numpoints,length)
        rot = self.convertdegtorad(init_rot)
        points = self.transformpoints(points,rot)
        data = self.simulatebuoyancy(points, init_velocity, init_rot)
        self.writedata(data)

    def simulatebuoyancy(self, points, init_velocity, init_rot):
        time = 0.0
        # Velocity and Angular Velocity are vectors
        velocity = init_velocity[0]
        angular_velocity = init_velocity[1]
        densityratio = self._fluid_density/self._obj_density
        acceleration = Vector(0,0-9.8) #assuming obj spawns in air
        data = {}
        times_list = []
        forces_list = []
        torques_list = []
        points_list = []
        rotation = init_rot

        while time < self._time_end:
            torque_sum = Vector(0,0,0)
            force_sum = Vector(0,0,0)

            newpoints = []
            for point in points:
                force_gravity = Vector(0,0,-(self.obj_mass/float(len(points)))*9.8)
                force_buoyancy = Vector(0,0,(self.obj_mass/float(len(points)))*densityratio*9.8)
                # assume point underwater
                lin_drag_coe = self._water_linear_damping
                ang_drag_coe = self._water_angular_damping
                
                # Cannot change points in iteration, points are connected
                x = velocity.x*self._time_interval + 0.5*acceleration.x*pow(self._time_interval,2.0) + point.x
                y = velocity.y*self._time_interval + 0.5*acceleration.y*pow(self._time_interval,2.0) + point.y
                z = velocity.z*self._time_interval + 0.5*acceleration.z*pow(self._time_interval,2.0) + point.z

                #check if point above water, then set air medium parameters
                if point.z > self._water_plane_z:
                    force_buoyancy.z = 0.0
                    lin_drag_coe = self._air_linear_damping
                    ang_drag_coe = self._air_angular_damping

                # Calculate drag variables

                trandragx = lin_drag_coe*(pow(velocity.x,2.0)/float(len(points)))
                trandragy = lin_drag_coe*(pow(velocity.y,2.0)/float(len(points)))
                trandragz = lin_drag_coe*(pow(velocity.z,2.0)/float(len(points)))
                force_linear_drag = Vector(trandragx,trandragy,trandragz)*(-velocity/velocity.magnitude)
                displacement_vector = self.getdispvector(point,points[len(points)/2])
                v_tan = angular_velocity.cross(displacement_vector)
                fad_x = (ang_drag_coe*pow(v_tan.x,2.0))/float(len(points))
                fad_y = (ang_drag_coe*pow(v_tan.y,2.0))/float(len(points))
                fad_z = (ang_drag_coe*pow(v_tan.z,2.0))/float(len(points))
                force_angular_drag = Vector(fad_x,fad_y,fad_z)*(-v_tan/v_tan.magnitude) 
                # Calculate total force and torque on point
                force_net = force_buoyancy + force_gravity + force_linear_drag 
                torque = displacement_vector.cross(force_net) 
                torque_drag = displacement_vector.cross(force_angular_drag)
                torque_net = torque + torque_drag
                
                # Summation: End of iterations will give the net torque and force on the entire object
                torque_sum += torque_net
                force_sum += force_net
                newpoint = Vector(x,y,z)
                newpoints.append(newpoint)

            acceleration = force_sum/self._rodmass
            angular_acceleration = torque_sum/self._moment_inertia
            velocity = self.updatelinearvelocity(velocity, self._time_interval, acceleration)
            angular_velocity = self.updateangularvelocity(angular_velocity, self._time_interval, angular_acceleration)
            d_rotation = self.updaterotation(angular_velocity, self._time_interval, angular_acceleration)
            newpoints = self.transformpoints(newpoints,d_rotation)
            # It is now safe to update the points to the new calculated set of points
            points = newpoints
            rotation = rotation + d_rotation

            times_list.append(time)
            points_list.append(points)
            torques_list.append(torque_sum)
            forces_list.append(force_sum)

            data = {
                'time':times_list,
                'points':points_list,
                'torque':torques_list,
                'force':forces_list
                }
            time += self._time_interval

            return data

    def writedata(self, data):
        with open('3DBuoyancy.csv', 'wb') as csvfile:
            writer = csv.DictWriter(csvfile, data.keys())
            writer.writerheader()
            times = data['time']
            points = data['points']
            forces = data['force']
            torques = data['torques']

            i = 0
            while i <len(times):
                row = {
                        'force':forces[i],
                        'torque':torques[i],
                        'location':points[i],
                        'time':times[i]
                        }
                writer.writerow(row)
                i+=1


    def getdispvector(self, point, origin):
        r_x = point.x - origin.x
        r_y = point.y - origin.y
        r_z = point.z - origin.z
        return Vector(r_x,r_y,r_y)

    def updatelinearvelocity(self, velocity, time, acceleration):
        vx = velocity.x + acceleration.x*time
        vy = velocity.y + acceleration.y*time
        vz = velocity.z + acceleration.z*time
        return Vector(vx,vy,vz)
    
    def updateangularvelocity(self, angular_velocity, time, angular_acceleration):
        wx = angular_velocity.x + angular_acceleration.x*time
        wy = angular_velocity.y + angular_acceleration.y*time
        wz = angular_velocity.z + angular_acceleration.z*time
        return Vector(wx,wy,wz)

    def updaterotation(self, angular_velocity, time, angular_acceleration):
        theta_x = (angular_velocity.x*time + (0.5*angular_acceleration.x*pow(time,2.0)))%(2*pi)
        theta_y = (angular_velocity.y*time + (0.5*angular_acceleration.y*pow(time,2.0)))%(2*pi)
        theta_z = (angular_velocity.z*time + (0.5*angular_acceleration.z*pow(time,2.0)))%(2*pi)
        return Vector(theta_x, theta_y, theta_z)

    def generatepoints(self, numpoints, length):
        points = []
        if numpoints < 3:
            print 'ERROR: invalid number of points.'
            sys.exit(0),'\n'
        if numpoints % 2 == 0:
            numpoints = numpoints + 1
            i = 0.0
            # Rod is being laid on the axis axis
            while i < numpoints:
                x = -rodLength/2.0 + ((rodLength*i)/(float(numpoints)-1.0))
                point = Vector(x,0,0)
                i+=1.0
        return points
    
    def convertdegtorad(self, radians):
        conversion = 57.2957795
        if isinstance(radians,Vector):
            radians.x = (radians.x * conversion)%(2*pi)
            radians.y = (radians.y * conversion)%(2*pi)
            radians.z = (radians.z * conversion)%(2*pi)
            return radians
        else: # assuming float
            return (radians * conversion)%(2*pi)

    def convertdegtorad(self, degrees):
        conversion = (1.0/57.2957795)
        if isinstance(degrees,Vector):
            degrees.x = (degrees.x * conversion)%(360)
            degrees.y = (degrees.y * conversion)%(360)
            degrees.z = (degrees.z * conversion)%(360)
            return degrees
        else: # assuming float
            return (degrees * conversion)%(360)
    
    def transformpoints(self, points, rotation):
        # Performs a 3D rotation counter clockwise in the x,y,z direction.
        rod_origin = points[len(points)/2]

        x_dist = obj_origin.x
        y_dist = obj_origin.y
        z_dist = obj_origin.z

        #move point to origin
        for point in points:
            x = point.x - x_dist
            y = point.y - y_dist
            z = point.z - z_dist
            v = Vector(x,y,z)
            v = self.xrotation(v,rotation.x)
            v = self.yrotation(v,rotation.y)
            v = self.zrotation(v,rotation.z)
            point = v

    def xrotation(self, vec, theta):
    # rotation around x axis:
        vec.x = vec.x
        vec.y = vec.y*cos(theta) + vec.z*sin(theta)
        vec.z = vec.z*cos(theta) - vec.y*sin(theta)
        return vec

    def yrotation(self, vec, theta):
        # rotation around y axis:
        vec.x = vec.x*cos(theta) - vec.z*sin(theta)
        vec.y = vec.y
        vec.z = vec.z*cos(theta) + vec.x*sin(theta)
        return vec

    def zrotation(self, vec, theta):
        # rotation around z axis:
        vec.x = vec.x*cos(theta) + vec.y*sin(theta)
        vec.y = vec.y*cos(theta) - vec.x*sin(theta)
        vec.z = vec.z
        return vec




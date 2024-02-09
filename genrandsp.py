import porespy as ps
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class generate_spherepack:
    def __init__(self, dimens, Rs, min_throat, overlapping=False, scale_factor=1, loop_limits=10000):
        self.dx = dimens[0]
        self.dy = dimens[1]
        self.dz = dimens[2]
        self.Rs = Rs
        self.rmax = np.max(self.Rs)
        self.rmin = np.min(self.Rs)
        self.min_throat = min_throat
        self.overlapping = overlapping
        self.scale_factor = scale_factor
        self.loop_limits = loop_limits
        self.loop_counter = 0
        self.sphere_counter = 0
        self.sphere_coords = []
        self.volume = self.dx*self.dy*self.dz
        
    def is_overlapping(self, sphere1, sphere2, min_throat):
        allowed_dist = sphere1[3]+sphere2[3]+min_throat
        distance = np.sqrt((sphere1[0]-sphere2[0])**2 + (sphere1[1]-sphere2[1])**2 + (sphere1[2]-sphere2[2])**2)
        if distance<allowed_dist:
            return True
        else: 
            return False
    
    def is_overlapping_all(self, x, y, z, r):
        # Duplicates generated sphere
        temp = np.tile([x,y,z,r], len(self.sphere_coords)).reshape([len(self.sphere_coords), 4])
        minths = np.tile(self.min_throat, len(self.sphere_coords)).reshape([len(self.sphere_coords), 1])
        # Compares generated sphere with all previously generated spheres to determine any overlapping
        overlapped = list(map(self.is_overlapping, self.sphere_coords, temp, minths)) # inputs ([100,4] , [100,4], min_throat)
                    

    def calc_porosity(self, sphere_coords, volume):
        sphere_coords = np.array(sphere_coords)
        all_radius = sphere_coords[:,3]
        return 1-np.sum(4/3*np.pi*(all_radius**3))/volume
    
    def generate_rand_sphere(self):
        r = np.random.choice(self.Rs)
        x = np.round(np.random.uniform(0+self.rmax,self.dx-self.rmax), decimals=6)
        y = np.round(np.random.uniform(0+self.rmax,self.dy-self.rmax), decimals=6)
        z = np.round(np.random.uniform(0+self.rmax,self.dz-self.rmax), decimals=6)
        return x, y ,z, r
    
    def update_loop_counter(self, reset=False):
        self.loop_counter += 1
        if reset==True:
            self.loop_counter = 0
        
    def update_sphere_counter(self, reset=False):
        self.sphere_counter += 1
        if reset==True:
            self.sphere_counter = 0
        
    def print_progress(self, timestep=500):
        if self.loop_counter%timestep==0:
            print(f"Loop #{self.loop_counter:0>5}, n_spheres={self.sphere_counter:0>5}, porosity={calc_porosity(self.sphere_coords, self.dx*self.dy*self.dz):.4f}")

    def reset_all(self):
        self.sphere_coords = []
        self.update_loop_counter(reset=True)
        self.update_sphere_counter(reset=True)

    # Generation methods
    def method_I(self, porosity):
        self.generate_rand_sphere()
        self.update_loop_counter()
        self.update_sphere_counter()
        
        while self.calc_porosity(self.sphere_coords, self.volume)>=porosity:
            if self.loop_counter==self.loop_limits:
                break
            
            if loop_counter%1000:
                r_range.pop()
            if overlapping==False:
                for r in list(reversed(r_range)):

                    self.is_overlapping_all(self.generate_rand_sphere())
                    # Duplicates generated sphere
                    temp = np.tile([x,y,z,r], len(sphere_coords)).reshape([len(sphere_coords), 4])
                    minths = np.tile(min_throat, len(sphere_coords)).reshape([len(sphere_coords), 1])
                    # Compares generated sphere with all previously generated spheres to determine any overlapping
                    overlapped = list(map(is_overlapping, sphere_coords, temp, minths)) # inputs ([100,4] , [100,4], min_throat)
                    # If there is no overlapping, add to the list
                    if all(not x for x in overlapped):
                        self.sphere_coords.append([x,y,z,r])
                        self.sphere_counter+=1
                        break



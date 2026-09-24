import porespy as ps
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def is_overlapping(sphere1, sphere2, min_throat):
    allowed_dist = sphere1[3]+sphere2[3]+min_throat
    distance = np.sqrt((sphere1[0]-sphere2[0])**2 + (sphere1[1]-sphere2[1])**2 + (sphere1[2]-sphere2[2])**2)
    if distance<allowed_dist:
        return True
    else: 
        return False


def calc_porosity(sphere_coords, volume):
    sphere_coords = np.array(sphere_coords)
    all_radius = sphere_coords[:,3]
    return 1-np.sum(4/3*np.pi*(all_radius**3))/volume
    

def generate_rand_spherepack_I(porosity, dx, dy, dz, rmin, rmax, min_throat, overlapping=False, scale_factor=1, loop_limits=10000):
    r"""
    Generate spheres until a certain porosity criteria is met

    output : array_like
        returns a numpy array of size (num_spheres, 4)
        each column represents x, y, z, radius of all spheres, respectively.
    
    parameters:
    ==================
    porosity : float
        desired porosity for the sphere pack, unfortunately it's not accurate yet
    dx, dy, dz : digits
        The length of the sphere pack in each direction.
    rmin, rmax : digits
        minimum and maximum radius of generated pores. 
        Set both parameters to the same value for same size spheres.
    min_throat : int or float
        minimum space between spheres. Doesn't matter if overlapping is True.
    overlapping : True/False
        Determines whether generated spheres overlap or not.
    scale_factor : int or float
        Scales the generated spheres for creating a 3D image.
        The values for dx,dy,dz and radius are relative. When converting to a 3D image, the radius of spheres needs to be an integer.
        Hence the min radius has to be at least one. The larger the scale factor, the higher the quality of sphere pack image, but also higher image size.
    

    """
    # Counters for spheres and number of loops
    sphere_counter = 1
    loop_counter = 0
    # define empty list
    sphere_coords = []
    # create one random sphere and add it to the list
    r = np.random.randint(rmin, rmax)
    x = np.round(np.random.uniform(0+r,dx-r), decimals=5)
    y = np.round(np.random.uniform(0+r,dy-r), decimals=5)
    z = np.round(np.random.uniform(0+r,dz-r), decimals=5)
    sphere_coords.append([x,y,z,r])
    # Start a loop to generate spheres
    while calc_porosity(sphere_coords, dx*dy*dz)>=porosity:
        loop_counter += 1
        # Breaks the loop
        if loop_counter==loop_limits:
            break
        # print current condition every 500 loops
        if loop_counter%500==0:
            print(f"Loop #{loop_counter:0>5}, n_spheres={sphere_counter:0>5}, poro={calc_porosity(sphere_coords, dx*dy*dz):.3f}")
        # create random location
        x = np.round(np.random.uniform(0+rmax,dx-rmax), decimals=5)
        y = np.round(np.random.uniform(0+rmax,dy-rmax), decimals=5)
        z = np.round(np.random.uniform(0+rmax,dz-rmax), decimals=5)
        r = np.random.randint(rmin, rmax)
        # It creates sphere with largest radius and reduces the radius until it doesn't overlap with any other sphere.
        r_range = np.arange(rmin, rmax+1)
        if loop_counter%1000:
            r_range.pop()
        if overlapping==False:
            for r in list(reversed(r_range)):
                # Duplicates generated sphere
                temp = np.tile([x,y,z,r], len(sphere_coords)).reshape([len(sphere_coords), 4])
                minths = np.tile(min_throat, len(sphere_coords)).reshape([len(sphere_coords), 1])
                # Compares generated sphere with all previously generated spheres to determine any overlapping
                overlapped = list(map(is_overlapping, sphere_coords, temp, minths)) # inputs ([100,4] , [100,4], min_throat)
                # If there is no overlapping, add to the list
                if all(not x for x in overlapped):
                    sphere_coords.append([x,y,z,r])
                    sphere_counter+=1
                    break
        # if overlapping==False:
        #     # Duplicates generated sphere
        #     temp = np.tile([x,y,z,r], len(sphere_coords)).reshape([len(sphere_coords), 4])
        #     minths = np.tile(min_throat, len(sphere_coords)).reshape([len(sphere_coords), 1])
        #     # Compares generated sphere with all previously generated spheres to determine any overlapping
        #     overlapped = list(map(is_overlapping, sphere_coords, temp, minths)) # inputs ([100,4] , [100,4], min_throat)
        #     # If there is no overlapping, add to the list
        #     if all(not x for x in overlapped):
        #         sphere_coords.append([x,y,z,r])
        #         sphere_counter+=1
        # If overlapping is True/allowed then it just creates a random sphere and adds it to the list.
        elif overlapping==True:
              r = np.random.randint(rmin, rmax)
              sphere_coords.append([x,y,z,r])
              sphere_counter+=1
    print(f"Loop #{loop_counter:0>5}, n_spheres={sphere_counter:0>5}, poro={calc_porosity(sphere_coords, dx*dy*dz):.3f}")
    return np.array(sphere_coords)


def generate_rand_spherepack_II(num_spheres, dx, dy, dz, rmin, rmax, min_throat, overlapping=False, scale_factor=1, loop_limits=20000):
    r"""
    Generate certain number of spheres.

    output : array_like
        returns a numpy array of size (num_spheres, 4)
        each column represents x, y, z, radius of all spheres, respectively.
    
    parameters:
    ==================
    num_spheres : int
        number of spheres to be generated
    dx, dy, dz : digits
        The length of the sphere pack in each direction.
    rmin, rmax : digits
        minimum and maximum radius of generated pores. 
        Set both parameters to the same value for same size spheres.
    min_throat : int or float
        minimum space between spheres. Doesn't matter if overlapping is True.
    overlapping : True/False
        Determines whether generated spheres overlap or not.
    scale_factor : int or float
        Scales the generated spheres for creating a 3D image.
        The values for dx,dy,dz and radius are relative. 
        When converting to a 3D image, the radius of spheres needs to be an integer. 
        Hence the min radius has to be at least one. 
        The larger the scale factor, the higher the quality of sphere pack image, but also higher image size.
    

    """
    # Counters for spheres and number of loops
    sphere_counter = 1
    loop_counter = 0
    # define empty list
    sphere_coords = []
    # create one random sphere and add it to the list
    r = np.random.randint(rmin, rmax)
    x = np.round(np.random.uniform(0+r,dx-r), decimals=5)
    y = np.round(np.random.uniform(0+r,dy-r), decimals=5)
    z = np.round(np.random.uniform(0+r,dz-r), decimals=5)
    sphere_coords.append([x,y,z,r])
    # Start a loop to generate spheres
    while sphere_counter<=num_spheres:
        loop_counter += 1
        # Breaks the loop
        if loop_counter==loop_limits:
            break
        # print current condition every 500 loops
        if loop_counter%1000==0:
            print(f"Loop #{loop_counter:0>5}, n_spheres={sphere_counter:0>5}, poro={calc_porosity(sphere_coords, dx*dy*dz):.3f}")
        # create random location
        x = np.round(np.random.uniform(0+rmax,dx-rmax), decimals=5)
        y = np.round(np.random.uniform(0+rmax,dy-rmax), decimals=5)
        z = np.round(np.random.uniform(0+rmax,dz-rmax), decimals=5)
        # It creates sphere with largest radius and reduces the radius until it doesn't overlap with any other sphere.
        r_range = np.arange(rmin, rmax+1)
        if overlapping==False:
            for r in list(reversed(r_range)):
                # Duplicates generated sphere
                temp = np.tile([x,y,z,r], len(sphere_coords)).reshape([len(sphere_coords), 4])
                minths = np.tile(min_throat, len(sphere_coords)).reshape([len(sphere_coords), 1])
                # Compares generated sphere with all previously generated spheres to determine any overlapping
                overlapped = list(map(is_overlapping, sphere_coords, temp, minths)) # inputs ([100,4] , [100,4], min_throat)
                # If there is no overlapping, add to the list
                if all(not x for x in overlapped):
                    sphere_coords.append([x,y,z,r])
                    sphere_counter+=1
                    break
        # If overlapping is True/allowed then it just creates a random sphere and adds it to the list.
        elif overlapping==True:
              r = np.random.randint(rmin, rmax)
              sphere_coords.append([x,y,z,r])
              sphere_counter+=1
    print(f"Loop #{loop_counter:0>5}, n_spheres={sphere_counter:0>5}, poro={calc_porosity(sphere_coords, dx*dy*dz):.3f}")
    return np.array(sphere_coords)



def generate_rand_spherepack_III(num_spheres, dx, dy, dz, Rs, min_throat, overlapping=False, scale_factor=1):
    r"""
    Generate spheres until certain numbers of spheres are generated,
    sphere radiuses have to be explicitly imported through the Rs parameter.

    output : array_like
        returns a numpy array of size (num_spheres, 4)
        each column represents x, y, z, radius of all spheres, respectively.
    
    parameters:
    ==================
    num_spheres : int
        number of spheres to be generated
    dx, dy, dz : digits
        The length of the sphere pack in each direction.
    Rs : array_like
        Contains radius of generated pores. 
        Set only one value to the array for same size spheres.
    min_throat : int or float
        minimum space between spheres. Doesn't matter if overlapping is True.
    overlapping : True/False
        Determines whether generated spheres overlap or not.
    scale_factor : int or float
        Scales the generated spheres for creating a 3D image.
        The values for dx,dy,dz and radius are relative. 
        When converting to a 3D image, the radius of spheres needs to be an integer. 
        Hence the min radius has to be at least one. 
        The larger the scale factor, the higher the quality of sphere pack image, but also higher image size.
    

    """
    # Counters for spheres and number of loops
    sphere_counter = 1
    loop_counter = 0
    # define empty list
    sphere_coords = []
    # create one random sphere and add it to the list
    r = np.random.choice(Rs)
    x = np.round(np.random.uniform(0+r,dx-r), decimals=5)
    y = np.round(np.random.uniform(0+r,dy-r), decimals=5)
    z = np.round(np.random.uniform(0+r,dz-r), decimals=5)
    sphere_coords.append([x,y,z,r])
    rmin, rmax = np.min(Rs), np.max(Rs)
    # Start a loop to generate spheres
    while num_spheres<=sphere_counter:
        loop_counter += 1
        # print current condition every 500 loops
        if loop_counter%500==0:
            print(f"Loop #{loop_counter:0>5}, n_spheres={sphere_counter:0>5}, poro={calc_porosity(sphere_coords, dx*dy*dz):.3f}")
        # create random location
        x = np.round(np.random.uniform(0+rmax,dx-rmax), decimals=5)
        y = np.round(np.random.uniform(0+rmax,dy-rmax), decimals=5)
        z = np.round(np.random.uniform(0+rmax,dz-rmax), decimals=5)
        # It creates sphere with largest radius and reduces the radius until it doesn't overlap with any other sphere.
        r = np.random.choice(Rs)
        if overlapping==False:
            # Duplicates generated sphere
            temp = np.tile([x,y,z,r], len(sphere_coords)).reshape([len(sphere_coords), 4])
            minths = np.tile(min_throat, len(sphere_coords)).reshape([len(sphere_coords), 1])
            # Compares generated sphere with all previously generated spheres to determine any overlapping
            overlapped = list(map(is_overlapping, sphere_coords, temp, minths)) # inputs ([100,4] , [100,4], min_throat)
            # If there is no overlapping, add to the list
            if all(not x for x in overlapped):
                sphere_coords.append([x,y,z,r])
                sphere_counter+=1
        # If overlapping is True/allowed then it just creates a random sphere and adds it to the list.
        elif overlapping==True:
              r = np.random.randint(rmin, rmax)
              sphere_coords.append([x,y,z,r])
              sphere_counter+=1
    print(f"Loop #{loop_counter:0>5}, n_spheres={sphere_counter:0>5}, poro={calc_porosity(sphere_coords, dx*dy*dz):.3f}")
    return np.array(sphere_coords)


def generate_rand_spherepack_IV(num_spheres, dx, dy, dz, Rs, min_throat, overlapping=False, scale_factor=1):
    r"""
    Generate spheres until certain numbers of spheres are generated,
    sphere radiuses have to be explicitly imported through the Rs parameter.

    output : array_like
        returns a numpy array of size (num_spheres, 4)
        each column represents x, y, z, radius of all spheres, respectively.
    
    parameters:
    ==================
    num_spheres : int
        number of spheres to be generated
    dx, dy, dz : digits
        The length of the sphere pack in each direction.
    Rs : array_like
        Contains radius of generated pores. 
        Set only one value to the array for same size spheres.
    min_throat : int or float
        minimum space between spheres. Doesn't matter if overlapping is True.
    overlapping : True/False
        Determines whether generated spheres overlap or not.
    scale_factor : int or float
        Scales the generated spheres for creating a 3D image.
        The values for dx,dy,dz and radius are relative. 
        When converting to a 3D image, the radius of spheres needs to be an integer. 
        Hence the min radius has to be at least one. 
        The larger the scale factor, the higher the quality of sphere pack image, but also higher image size.
    

    """
    # Counters for spheres and number of loops
    sphere_counter = 1
    loop_counter = 0
    # define empty list
    sphere_coords = []
    # create one random sphere and add it to the list
    r = np.random.choice(Rs)
    x = np.round(np.random.uniform(0+r,dx-r), decimals=5)
    y = np.round(np.random.uniform(0+r,dy-r), decimals=5)
    z = np.round(np.random.uniform(0+r,dz-r), decimals=5)
    sphere_coords.append([x,y,z,r])
    rmin, rmax = np.min(Rs), np.max(Rs)
    # Start a loop to generate spheres
    while num_spheres<=sphere_counter:
        loop_counter += 1
        # print current condition every 500 loops
        if loop_counter%500==0:
            print(f"Loop #{loop_counter:0>5}, n_spheres={sphere_counter:0>5}, poro={calc_porosity(sphere_coords, dx*dy*dz):.3f}")
        # create random location
        x = np.round(np.random.uniform(0+rmax,dx-rmax), decimals=5)
        y = np.round(np.random.uniform(0+rmax,dy-rmax), decimals=5)
        # z = np.round(np.random.uniform(0+rmax,dz-rmax), decimals=5)
        r = np.random.choice(Rs)
        z = dz-r
        if overlapping==False:
            # Duplicates generated sphere
            temp = np.tile([x,y,z,r], len(sphere_coords)).reshape([len(sphere_coords), 4])
            minths = np.tile(min_throat, len(sphere_coords)).reshape([len(sphere_coords), 1])
            # Compares generated sphere with all previously generated spheres to determine any overlapping
            overlapped = list(map(is_overlapping, sphere_coords, temp, minths)) # inputs ([100,4] , [100,4], min_throat)
            # If there is no overlapping, add to the list
            while z>0 and any(overlapped): # loop until overlapping==True and still in bounding box
                if all(not x for x in overlapped):
                    sphere_coords.append([x,y,z,r])
                    sphere_counter+=1
                else:
                    z -= dz//100       

        # If overlapping is True/allowed then it just creates a random sphere and adds it to the list.
        elif overlapping==True:
              r = np.random.randint(rmin, rmax)
              sphere_coords.append([x,y,z,r])
              sphere_counter+=1
    print(f"Loop #{loop_counter:0>5}, n_spheres={sphere_counter:0>5}, poro={calc_porosity(sphere_coords, dx*dy*dz):.3f}")
    return np.array(sphere_coords)
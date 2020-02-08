import os
import sys
import numpy as np
import math
import random as rand

import util_plot
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../lib_mag')
import MagTenseSource

class Tiles():
        def __init__(self, n):
                # Initialization of arrays for specific tile parameters
                # Input to Fortran derived type MagTile     
                self.center_pos = np.zeros(shape=(n,3), dtype=np.float64, order='F') # r0, theta0, z0
                self.dev_center = np.zeros(shape=(n,3), dtype=np.float64, order='F') # dr, dtheta, dz
                self.size = np.zeros(shape=(n,3), dtype=np.float64, order='F') # a, b, c
                self.vertices = np.zeros(shape=(n,3,4), dtype=np.float64, order='F') # v1, v2, v3, v4 as column vectors
                self.M = np.zeros(shape=(n,3), dtype=np.float64, order='F') # Mx, My, Mz
                self.u_ea = np.zeros(shape=(n,3), dtype=np.float64, order='F') # Easy axis
                self.u_oa1 = np.zeros(shape=(n,3), dtype=np.float64, order='F')
                self.u_oa2 = np.zeros(shape=(n,3), dtype=np.float64, order='F')
                self.mu_r_ea = np.ones(shape=(n), dtype=np.float64, order='F')
                self.mu_r_oa = np.ones(shape=(n), dtype=np.float64, order='F')
                self.M_rem = np.zeros(shape=(n), dtype=np.float64, order='F')
                self.tile_type = np.ones(n, dtype=np.int32, order='F') # 1 = cylinder, 2 = prism, 3 = circ_piece, 4 = circ_piece_inv, 5 = tetrahedron, 10 = ellipsoid
                self.offset = np.zeros(shape=(n,3), dtype=np.float64, order='F') # offset of global coordinates
                self.rot = np.zeros(shape=(n,3), dtype=np.float64, order='F')
                self.color = np.zeros(shape=(n,3), dtype=np.float64, order='F')
                self.magnetic_type = np.ones(n, dtype=np.int32, order='F') # 1 = hard magnet, 2 = soft magnet
                self.stfcn_index = np.ones(shape=(n), dtype=np.int32, order='F') # default index into the state function
                self.incl_it = np.ones(shape=(n), dtype=np.int32, order='F') # if equal to zero the tile is not included in the iteration
                self.use_sym = np.zeros(shape=(n), dtype=np.int32, order='F') # whether to exploit symmetry
                self.sym_op = np.ones(shape=(n,3), dtype=np.float64, order='F') # 1 for symmetry and -1 for anti-symmetry respectively to the planes
                self.M_rel = np.zeros(shape=(n), dtype=np.float64, order='F')                
                
                self.grid_pos = np.zeros(shape=(n,3), dtype=np.float64, order='F') # positions in the grid
                self.n = n

        def __str__(self):
                result = ("{n} tiles are present in the setup:\n").format(n = self.n)
                for i in range(self.n):
                        result = result + ("Tile {i} at grid position ({grid_x},{grid_y},{grid_z}) with coordinates x={x}, y={y}, z={z}.\n").format(\
                        i = i, grid_x = self.grid_pos[i][0], grid_y = self.grid_pos[i][1], grid_z=self.grid_pos[i][2], \
                        x=self.offset[i][0], y = self.offset[i][1], z=self.offset[i][2])
                return result

        def get_n(self):
                return self.n

        def set_grid_pos(self, grid_pos):
                for i,pos in enumerate(grid_pos):
                        self.set_grid_pos_i(pos,i)

        def set_grid_pos_i(self, grid_pos, i):
                self.grid_pos[i] = grid_pos
             
        def set_size(self, sizes):
                if isinstance(sizes[0], int) or isinstance(sizes[0], float):
                        self.size[:] = sizes
                else:
                        for i,size in enumerate(sizes):
                                self.set_size_i(size,i)

        def set_size_i(self, size, i):
                self.size[i] = size
        
        def get_size(self, i):
                return self.size[i]
        
        def set_center_pos(self, center_positions):
                for i,center_pos in enumerate(center_positions):
                        self.set_center_pos_i(center_pos,i)

        def set_center_pos_i(self, center_pos, i):
                self.center_pos[i] = center_pos
        
        def get_center_pos(self, i):
                return self.center_pos[i]
        
        def set_dev_center(self, devs):
                if isinstance(devs[0], int) or isinstance(devs[0], float):
                        self.dev_center[:] = devs
                else:
                        for i,dev_center in enumerate(devs):
                                self.set_dev_center_i(dev_center,i)

        def set_dev_center_i(self, dev_center, i):
                self.dev_center[i] = dev_center
        
        def get_dev_center(self, i):
                return self.dev_center[i]

        def set_vertices(self, n_vertices):
                if n_vertices.shape == (3,4):
                        self.vertices[:] = n_vertices
                elif n_vertices.shape == (4,3):
                        self.vertices[:] = n_vertices.transpose()
                else:
                        for i,vertices in enumerate(n_vertices):
                                self.set_vertices_i(vertices,i)

        def set_vertices_i(self, vertices, i):
                if vertices.shape == (3,4):
                        self.vertices[i] = vertices
                elif vertices.shape == (4,3):
                        self.vertices[i] = vertices.transpose()
                else:
                        print("Four 3-dimensional vertices have to be defined!")
                
        
        def get_vertices(self, i):
                return self.vertices[i]
        
        def set_tile_type(self, tile_types):
                if isinstance(tile_types, int) or isinstance(tile_types, float):
                        self.tile_type[:] = tile_types
                else:
                        for i,tile_type in enumerate(tile_types):
                                self.set_tile_type_i(tile_type,i)

        def set_tile_type_i(self, tile_type, i):
                self.tile_type[i] = tile_type
        
        def get_tile_type(self, i):
                return self.tile_type[i]

        def set_offset_i(self, offset, i):
                self.offset[i] = offset

        def get_offset(self, i):
                return self.offset[i]
        
        def set_rotation_i(self, rotation, i):
                self.rot[i] = rotation

        def get_rotation(self, i):
                return self.rot[i]
        
        def set_easy_axis(self, easy_axis):
                for i,ea in enumerate(easy_axis):
                        self.u_ea[i] = np.around(ea, decimals=9)
                        self.M[i] = self.M_rem[i] * self.u_ea[i]
                        oa_1 = np.array([easy_axis[i][1], -easy_axis[i][0], 0])
                        oa_1 = oa_1 / np.linalg.norm(oa_1)
                        self.u_oa1[i] = np.around(oa_1, decimals=9)
                        oa_2 = np.cross(self.u_ea[i], self.u_oa1)
                        self.u_oa2[i] = np.around(oa_2, decimals=9)

        def set_easy_axis_i(self, easy_axis, i):
                self.u_ea[i] = np.around(easy_axis, decimals=9)
                self.M[i] = self.M_rem[i] * self.u_ea[i]
        
        def set_oa1_i(self, other_axis, i):
                self.u_oa1[i] = np.around(other_axis, decimals=9)
        
        def set_oa2_i(self, other_axis, i):
                self.u_oa2[i] = np.around(other_axis, decimals=9)
        
        def set_mu_r_ea(self, mu):
                if isinstance(mu, int) or isinstance(mu, float):
                        self.mu_r_ea[:] = mu
                else:
                        for i,mu_i in enumerate(mu):
                                self.set_mu_r_ea_i(mu_i,i)
        
        def set_mu_r_ea_i(self, mu, i):
                self.mu_r_ea[i] = mu
        
        def set_mu_r_oa(self, mu):
                if isinstance(mu, int) or isinstance(mu, float):
                        self.mu_r_oa[:] = mu
                else:
                        for i,mu_i in enumerate(mu):
                                self.mu_r_oa(mu_i,i)

        def set_mu_r_oa_i(self, mu, i):
                self.mu_r_oa[i] = mu
        
        def set_remanence(self, M_rem):
                if isinstance(M_rem, int) or isinstance(M_rem, float):
                        self.M_rem[:] = M_rem
                else:
                        for i,M_rem_i in enumerate(M_rem):
                                self.set_remanence_i(M_rem_i,i)
        
        def set_remanence_i(self, M_rem, i):
                self.M_rem[i] = M_rem
        
        def set_mag_angle(self, mag_angles):
                for i,mag_angle in enumerate(mag_angles):
                        self.set_mag_angle_i(mag_angle,i)

        def set_mag_angle_rand(self):
                for i in range(self.n):
                        self.set_mag_angle_i([math.pi * rand.random(), 2*math.pi * rand.random()], i)                      

        def set_mag_angle_i(self, spherical_angles, i):
                # polar angle [0, pi], azimuth [0, 2*pi]
                if isinstance(spherical_angles, int) or isinstance(spherical_angles, float):
                        print("Azimuth and polar angle have to be set!\nExiting!")
                        exit()
                else:
                        polar_angle = spherical_angles[0]
                        azimuth = spherical_angles[1]
                        self.set_easy_axis_i([math.sin(polar_angle) * math.cos(azimuth), math.sin(polar_angle) * math.sin(azimuth), math.cos(polar_angle)], i)
                        self.set_oa1_i([math.sin(polar_angle) * math.sin(azimuth), math.sin(polar_angle) * (-math.cos(azimuth)), 0], i)
                        self.set_oa2_i([0.5*math.sin(2*polar_angle) * math.cos(azimuth), 0.5*math.sin(2*polar_angle) * math.sin(azimuth), -math.pow(math.sin(polar_angle),2)], i)

        def set_M(self, M, i):
                self.M[i] = M
        
        def get_M(self, i):
                return self.M[i]

        def set_color(self, color):
                if isinstance(color[0], int) or isinstance(color[0], float):
                        self.color[:] = color
                else:
                        for i,color_i in enumerate(color):
                                self.set_color_i(color_i,i)
        
        def set_color_i(self, color, i):
                self.color[i] = color

        def get_color(self, i):
                return self.color[i]

        def set_mag_type_i(self, mag_type, i):
                self.magnetic_type[i] = mag_type
  
        def set_as_Fe(self, idx):
                if isinstance(idx, int) or isinstance(idx, float):
                        self.set_as_Fe_i(idx)
                else:
                        for i in idx:
                                self.set_as_Fe_i(i)

        def set_as_Fe_i(self, i):
                self.set_remanence_i(0,i) # No remanence for iron
                self.set_mag_angle_i([math.pi/2, 0], i) # Default: Easy axis in direction of x-axis
                self.set_mu_r_ea_i(4000,i) # High relative permeability for iron
                self.set_mu_r_oa_i(4000,i) # High relative permeability for iron
                self.set_color_i([0,0.5,0.2],i) # Dark Green
                self.set_mag_type_i(2,i)
                self.set_M([1000,0,0],i)
        
        def add_tiles(self, n):
                self.center_pos = np.append(self.center_pos, np.zeros(shape=(n,3), dtype=np.float64, order='F'), axis = 0) # r0, theta0, z0
                self.dev_center = np.append(self.dev_center, np.zeros(shape=(n,3), dtype=np.float64, order='F'), axis = 0) # dr, dtheta, dz
                self.size = np.append(self.size, np.zeros(shape=(n,3), dtype=np.float64, order='F'), axis = 0) # a, b, c
                self.vertices = np.append(self.vertices, np.zeros(shape=(n,3,4), dtype=np.float64, order='F'), axis = 0) # v1, v2, v3, v4 as column vectors
                self.M = np.append(self.M, np.zeros(shape=(n,3), dtype=np.float64, order='F'), axis = 0) # Mx, My, Mz
                self.u_ea = np.append(self.u_ea, np.zeros(shape=(n,3), dtype=np.float64, order='F'), axis = 0) # Easy axis
                self.u_oa1 = np.append(self.u_oa1, np.zeros(shape=(n,3), dtype=np.float64, order='F'), axis = 0)
                self.u_oa2 = np.append(self.u_oa2, np.zeros(shape=(n,3), dtype=np.float64, order='F'), axis = 0)
                self.mu_r_ea = np.append(self.mu_r_ea, np.ones(shape=(n), dtype=np.float64, order='F'), axis = 0)
                self.mu_r_oa = np.append(self.mu_r_oa, np.ones(shape=(n), dtype=np.float64, order='F'), axis = 0)
                self.M_rem = np.append(self.M_rem, np.zeros(shape=(n), dtype=np.float64, order='F'), axis = 0)
                self.tile_type = np.append(self.tile_type, np.ones(n, dtype=np.int32, order='F'), axis = 0) # 1 = cylinder, 2 = prism, 3 = circ_piece, 4 = circ_piece_inv, 5 = tetrahedron, 10 = ellipsoid
                self.offset = np.append(self.offset, np.zeros(shape=(n,3), dtype=np.float64, order='F'), axis = 0) # offset of global coordinates
                self.rot = np.append(self.rot, np.zeros(shape=(n,3), dtype=np.float64, order='F'), axis = 0)
                self.color = np.append(self.color, np.zeros(shape=(n,3), dtype=np.float64, order='F'), axis = 0)
                self.magnetic_type = np.append(self.magnetic_type, np.ones(n, dtype=np.int32, order='F'), axis = 0) # 1 = hard magnet, 2 = soft magnet
                self.stfcn_index = np.append(self.stfcn_index, np.ones(shape=(n), dtype=np.int32, order='F'), axis = 0) # default index into the state function
                self.incl_it = np.append(self.incl_it, np.ones(shape=(n), dtype=np.int32, order='F'), axis = 0) # if equal to zero the tile is not included in the iteration
                self.use_sym = np.append(self.use_sym, np.zeros(shape=(n), dtype=np.int32, order='F'), axis = 0) # whether to exploit symmetry
                self.sym_op = np.append(self.sym_op, np.ones(shape=(n,3), dtype=np.float64, order='F'), axis = 0) # 1 for symmetry and -1 for anti-symmetry respectively to the planes
                self.M_rel = np.append(self.M_rel, np.zeros(shape=(n), dtype=np.float64, order='F'), axis = 0)                             
                self.grid_pos = np.append(self.grid_pos, np.zeros(shape=(n,3), dtype=np.float64, order='F'), axis = 0) # positions in the grid
                self.n = self.n + n

        def refinement_prism(self, idx):
                if isinstance(idx, int) or isinstance(idx, float):
                        self.refinement_prism_i(idx)
                else:
                        for i in idx:
                                self.refinement_prism_i(i)
        
        def refinement_prism_i(self, i):
                n = self.n
                self.add_tiles(7)
                center_point = self.offset[i]
                x_off, y_off, z_off = self.size[i]/4

                self.size[n:n+7] = self.size[i]/2
                self.size[i] = self.size[i]/2
                self.M[n:n+7] = self.M[i]
                self.M_rel[n:n+7] = self.M_rel[i]
                self.color[n:n+7] = self.color[i]
                self.magnetic_type[n:n+7] = self.magnetic_type[i]
                self.mu_r_ea[n:n+7] = self.mu_r_ea[i]
                self.mu_r_oa[n:n+7] = self.mu_r_oa[i]
                self.rot[n:n+7] = self.rot[i]
                self.tile_type[n:n+7] = self.tile_type[i]
                self.u_ea[n:n+7] = self.u_ea[i]
                self.u_oa1[n:n+7] = self.u_oa1[i]
                self.u_oa2[n:n+7] = self.u_oa2[i]

                self.offset[n] = center_point + [x_off, y_off, z_off]
                self.offset[n+1] = center_point + [-x_off, y_off, z_off]
                self.offset[n+2] = center_point + [x_off, -y_off, z_off]
                self.offset[n+3] = center_point + [-x_off, -y_off, z_off]
                self.offset[n+4] = center_point + [x_off, y_off, -z_off]
                self.offset[n+5] = center_point + [-x_off, y_off, -z_off]
                self.offset[n+6] = center_point + [x_off, -y_off, -z_off]
                self.offset[i] = center_point + [-x_off, -y_off, -z_off]
                
                

class Grid():
        def __init__(self, places, area):
                self.places = np.asarray(places)
                self.area = np.asarray(area)
                self.size_tile = self.area/self.places
                self.grid = np.zeros(self.places.tolist())

        def get_tiles(self, grid_positions = [], n_tiles = None):
                tiles = []
                # Fill grid randomly
                if not grid_positions:
                        if n_tiles is None:
                                n_tiles = rand.randrange(np.prod(self.places))
                        else:
                                n_tiles = n_tiles
                        grid_positions = []
                        
                        while n_tiles > 0:
                                new_pos = [rand.randrange(self.places[0]), rand.randrange(self.places[1]), rand.randrange(self.places[2])]
                                if new_pos in grid_positions:
                                        continue
                                else:
                                        grid_positions.append(new_pos)
                                        n_tiles = n_tiles - 1
                
                if len(grid_positions) == 0:
                        tiles = None
                else:
                        tiles = Tiles(len(grid_positions))
                        # Set grid position of tile
                        tiles.set_grid_pos(grid_positions)
                        # Set size of tile
                        tiles.set_size(self.size_tile)
                        # Set tile type: 2 = prism
                        tiles.set_tile_type(2)

                        for i,pos in enumerate(grid_positions):
                                if True in np.greater_equal(np.asarray(pos), self.places):
                                        print(("Desired position {} is not in the grid!").format(pos))
                                        exit()
                                self.grid[pos[0],pos[1],pos[2]] = 1

                                # Extract Cartesian coordinates of tile
                                tiles.set_offset_i(np.around((pos * self.size_tile) + self.size_tile/2, decimals=9),i)

                return tiles
        
        def set_eval_points(self, n_points, mode):
                counter = 0
                points = np.zeros(shape=(n_points[0]*n_points[1]*n_points[2],3), dtype=np.float64, order='F')
                if mode == "uniform":
                        seg = self.area/np.asarray(n_points)
                        for i in range(0,n_points[0]):
                                for j in range(0, n_points[1]):
                                        for k in range(0, n_points[2]):
                                                points[counter] = [i*seg[0]+seg[0]/2, j*seg[1]+seg[1]/2, k*seg[2]+seg[2]/2]
                                                counter = counter + 1
                elif mode == "center":
                        center = self.area/2
                        seg = (self.area[0]/10)/n_points[0]
                        seg_angle = math.pi/n_points[1]
                        seg_layer = self.area[2]/n_points[2]
                        for i in range(0,n_points[0]):
                                for j in range(0, n_points[1]):
                                        for k in range(0, n_points[2]):
                                                points[counter] = [i*seg, j*seg_angle, k*seg_layer] + center
                                                counter = counter + 1
                elif mode == "costumized":
                        pass
                else:
                        print("Please specify a valid area of interest!")
                        exit()
                return points

        def clear_grid(self):
                self.grid = np.zeros(self.places.tolist())

def get_average_magnetic_flux(H):
        norm = get_norm_magnetic_flux(H)
        average = sum(norm)/len(norm)
        return average

def get_p2p(H):
        norm = get_norm_magnetic_flux(H)
        return max(norm)-min(norm)

def get_norm_magnetic_flux(H):
        norm = np.zeros(shape=len(H))
        mu0 = 4*math.pi*1e-7 # vacuum permeability
        norm = [np.linalg.norm(H_point)*mu0 for H_point in H]
        return norm

def get_norm_magnetic_field(H):
        return [np.linalg.norm(H_point) for H_point in H]


def setup(places, area, n_tiles=0, filled_positions=None, mag_angles=[], eval_points=[20, 20, 5], eval_mode="uniform", B_rem=1.2):
        # Check format of input parameters
        if len(places) != 3:
                print("Format of number of possible magnets in each axis is not correct!")
                exit
        elif len(area) != 3:
                print("Format of area is not correct!")
                exit
        elif len(eval_points) != 3:
                print("Format of eval points is not correct!")
                exit
        meshgrid = Grid(places, area)
        # Extract coordinates of evaluation points
        points = meshgrid.set_eval_points(eval_points, eval_mode)
        # Fill grid with magnetic tiles
        if filled_positions is None:
                tiles = meshgrid.get_tiles(n_tiles=n_tiles)
        else:
                tiles = meshgrid.get_tiles(grid_positions=filled_positions)
        # Assign magnetization angles for tiles
        if tiles is not None:
                if not mag_angles:
                        # polar angle [0, pi], azimuth [0, 2*pi]
                        for _ in range(tiles.get_n()):
                                mag_angles.append([math.pi * rand.random(), 2*math.pi * rand.random()])

                for i in range(tiles.get_n()):        
                        if not mag_angles[i]:
                                mag_angles[i] = [math.pi * rand.random(), 2*math.pi * rand.random()]

                tiles.set_remanence(B_rem / (4*math.pi*1e-7))
                tiles.set_mag_angle(mag_angles)
                # Setting display color of magnets: red - [1, 0, 0] 
                tiles.set_color([1, 0, 0])     
        return tiles, points, meshgrid

# Function for running MagTense with the Fortran source code as Python module
def run_simulation(tiles, points, grid=None, plot=False, max_error=0.00001, max_it=500, iterate_solution=True, return_field=True, T = 300.):
        max_error = max_error # max relative error
        max_it = max_it # max number of iterations 
        iterate_solution = iterate_solution
        return_field = return_field
        T = T # temperature for the state function of iron (arbitrary here as we ignore temperature variation in the iron)
        data_stateFcn = np.genfromtxt(os.path.dirname(os.path.abspath(__file__)) + '/../util/data_stateFcn_saturated.csv', delimiter=';', dtype=np.float64)

        H, M_out, Mrel_out = \
                MagTenseSource.fortrantopythonio.runsimulation( centerpos=tiles.center_pos, dev_center=tiles.dev_center, \
                rect_size=tiles.size, vertices=tiles.vertices, mag=tiles.M, u_ea=tiles.u_ea, u_oa1=tiles.u_oa1, u_oa2=tiles.u_oa2, mu_r_ea=tiles.mu_r_ea, \
                mu_r_oa=tiles.mu_r_oa, mrem=tiles.M_rem, tiletype=tiles.tile_type, offset=tiles.offset, rotangles=tiles.rot, \
                color=tiles.color, magnettype=tiles.magnetic_type, statefunctionindex=tiles.stfcn_index, includeiniteration=tiles.incl_it, \
                exploitsymmetry=tiles.use_sym, symmetryops=tiles.sym_op, mrel=tiles.M_rel, pts=points, data_statefcn=data_stateFcn, \
                n_statefcn=1, t=T, maxerr=max_error, nitemax=max_it, iteratesolution=iterate_solution, returnsolution=return_field  )

        updated_tiles = tiles
        
        if iterate_solution is True:
                updated_tiles.M = M_out
                # updated_tiles.u_ea = u_ea_out
                # updated_tiles.u_oa1 = u_oa1_out
                # updated_tiles.u_oa2 = u_oa2_out
                updated_tiles.M_rel = Mrel_out
        
        if return_field is True:
                solution = H
        else:
                solution = None

        if plot==True:
                if grid is None:
                        util_plot.create_plot(updated_tiles, points, H)
                else:
                        util_plot.create_plot(updated_tiles, points, H, grid=grid)
        
        return updated_tiles, solution

def iterate_magnetization(tiles, max_error=0.00001, max_it=500, T=300.):
        max_error = max_error # max relative error
        max_it = max_it # max number of iterations
        T = T # temperature for the state function of iron (arbitrary here as we ignore temperature variation in the iron)
        data_stateFcn = np.genfromtxt(os.path.dirname(os.path.abspath(__file__)) + '/../util/data_stateFcn.csv', delimiter=';', dtype=np.float64)

        M_out, Mrel_out = MagTenseSource.fortrantopythonio.iteratetiles( centerpos=tiles.center_pos, dev_center=tiles.dev_center, \
                rect_size=tiles.size, vertices=tiles.vertices, mag=tiles.M, u_ea=tiles.u_ea, u_oa1=tiles.u_oa1, u_oa2=tiles.u_oa2, mu_r_ea=tiles.mu_r_ea, \
                mu_r_oa=tiles.mu_r_oa, mrem=tiles.M_rem, tiletype=tiles.tile_type, offset=tiles.offset, rotangles=tiles.rot, color=tiles.color, \
                magnettype=tiles.magnetic_type, statefunctionindex=tiles.stfcn_index, includeiniteration=tiles.incl_it, exploitsymmetry=tiles.use_sym, \
                symmetryops=tiles.sym_op, mrel=tiles.M_rel, data_statefcn=data_stateFcn, n_statefcn=1, t=T, maxerr=max_error, nitemax=max_it )

        updated_tiles = tiles
        updated_tiles.M = M_out
        updated_tiles.M_rel = Mrel_out

        return updated_tiles

def get_N_tensor(tiles, points):
        N = MagTenseSource.fortrantopythonio.getnfromtiles( centerpos=tiles.center_pos, dev_center=tiles.dev_center, \
                rect_size=tiles.size, vertices=tiles.vertices, mag=tiles.M, u_ea=tiles.u_ea, u_oa1=tiles.u_oa1, u_oa2=tiles.u_oa2, mu_r_ea=tiles.mu_r_ea, \
                mu_r_oa=tiles.mu_r_oa, mrem=tiles.M_rem, tiletype=tiles.tile_type, offset=tiles.offset, rotangles=tiles.rot, color=tiles.color, \
                magnettype=tiles.magnetic_type, statefunctionindex=tiles.stfcn_index, includeiniteration=tiles.incl_it, exploitsymmetry=tiles.use_sym, \
                symmetryops=tiles.sym_op, mrel=tiles.M_rel, pts=points )                   

        return N

def get_H_field(tiles, points, N=None):
        if N is None:
                N = np.zeros(shape=(tiles.get_n(),len(points),3,3), dtype=np.float64, order='F')
                useN = False
        else:
                useN = True

        H = MagTenseSource.fortrantopythonio.gethfromtiles( centerpos=tiles.center_pos, dev_center=tiles.dev_center, \
                rect_size=tiles.size, vertices=tiles.vertices, mag=tiles.M, u_ea=tiles.u_ea, u_oa1=tiles.u_oa1, u_oa2=tiles.u_oa2, mu_r_ea=tiles.mu_r_ea, \
                mu_r_oa=tiles.mu_r_oa, mrem=tiles.M_rem, tiletype=tiles.tile_type, offset=tiles.offset, rotangles=tiles.rot, color=tiles.color, \
                magnettype=tiles.magnetic_type, statefunctionindex=tiles.stfcn_index, includeiniteration=tiles.incl_it, exploitsymmetry=tiles.use_sym, \
                symmetryops=tiles.sym_op, mrel=tiles.M_rel, pts=points, n=N, useprovidedn=useN )                   

        return H

import os
import sys
import numpy as np
import math
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../source')
import MagTense
import MagTenseStandalone
import util_plot

def main():
    height = 0.01
    # Defining grid
    eval_offset = [0.45, 0.55, 0.025]
    places = [10, 10, 5]
    area = [1, 1, height*5]
    # Defining occupied places in grid
    filled_positions = [[4, 5, 2]]
    # Defining angle of magnetization in spherical coordinates (azimuth, polar angle) for each tile
    mag_angles = [[math.pi/2, math.pi/2]]

    # Optional parameters for setup: n_magnets, filled_positions, mag_angles, eval_points, eval_mode, B_rem
    (tiles, points, grid) = MagTense.setup(places, area, filled_positions=filled_positions, mag_angles=mag_angles)

    tiles.set_mu_r_ea(1.06)
    tiles.set_mu_r_oa(1.17)

    #Set costum points for comparison with COMSOL
    x_values = np.linspace(0,area[0],200)
    y_values = np.linspace(0,area[1],200)
    z_values = np.linspace(0,area[2],100)

    # tiles.set_as_Fe([1,2])
    # tiles.refinement_prism(1)
    # tiles.refinement_prism(10)
    
    struc_x = np.ones(len(x_values))
    pts_x = np.c_[x_values, struc_x*eval_offset[1], struc_x*eval_offset[2]]

    struc_y = np.ones(len(y_values))
    pts_y = np.c_[struc_y*eval_offset[0], y_values, struc_y*eval_offset[2]]

    struc_z = np.ones(len(z_values))
    pts_z = np.c_[struc_z*eval_offset[0], struc_z*eval_offset[1], z_values]

    # Standard parameters in settings: max_error=0.00001, max_it=500
    # All steps in one command - iterate_solution = True, return_field = True
    # N is not reused for calculating the H field
    (updated_tiles, H_x) = MagTense.run_simulation(tiles, pts_x, grid=grid, plot=False)
    (updated_tiles, H_y) = MagTense.run_simulation(tiles, pts_y, grid=grid, plot=False)
    (updated_tiles, H_z) = MagTense.run_simulation(tiles, pts_z, grid=grid, plot=False)

    #Plotting B-field for each axis through center
    for (pts,H_field) in [(x_values,H_x),(y_values,H_y),(z_values,H_z)]:
        fig = plt.figure()
        ax = fig.gca()
        ax.plot(pts, MagTense.get_norm_magnetic_flux(H_field), 'r*', label='MagTense')
        ax.legend()
        ax.set_xlabel('coordinate')
        ax.set_ylabel('B')
        plt.title("Validating magnetic field")
        plt.show()

main()

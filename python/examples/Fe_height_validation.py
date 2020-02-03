import os
import sys
import numpy as np
import math
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../source')
import MagTense
import MagTenseStandalone
import util_plot
import util_eval

def main():
    # Load reference points from COMSOL calculation
    COMSOL_eval_path = os.path.dirname(os.path.abspath(__file__)) + '/../util/evaluations/'
    # Defining grid
    places = [10, 10, 1]
    area = [1, 1, 0.01]
    # Defining occupied places in grid
    filled_positions = [[3, 3, 0], [3, 5, 0], [5, 3, 0], [5, 5, 0]]
    # Defining angle of magnetization in spherical coordinates (azimuth, polar angle) for each tile
    mag_angles = [[math.pi/2, math.pi/2], [], [], [math.pi/2, 2*math.pi]]

    # Optional parameters for setup: n_magnets, filled_positions, mag_angles, eval_points, eval_mode, B_rem
    (tiles, points, grid) = MagTense.setup(places, area, filled_positions=filled_positions, mag_angles=mag_angles)

    eval_offset = [0.55, 0.35, 0]
    (eval_points_z, B_h_COMSOL) = util_eval.load_COMSOL_eval('Stefan_iron_B_h.txt', eval_offset, COMSOL_eval_path)
    id_run = [8,739,1477,2169,2869,3571,4285,4979,5692,6366]

    tiles.set_mu_r_ea(1.06)
    tiles.set_mu_r_oa(1.17)

    tiles.set_as_Fe([1,2])

    B_h_MagTense = []
    # z_coor = np.linspace(-1,1,501)
    # struc = np.ones(len(z_coor))
    # pts = np.c_[struc*eval_offset[0], struc*eval_offset[1], z_coor]

    for z in range(1, 11):
        tiles.set_size([0.1, 0.1, z/100])
        z_coor = eval_points_z[(id_run[z-1]):id_run[z],2]
        struc = np.ones(len(z_coor))
        pts = np.c_[struc*eval_offset[0], struc*eval_offset[1], z_coor]
        # Standard parameters in settings: max_error=0.00001, max_it=500
        # All steps in one command - iterate_solution = True, return_field = True
        # N is not reused for calculating the H field
        (updated_tiles, H) = MagTense.run_simulation(tiles, pts, grid=grid, plot=False)
        #B_h_MagTense.append(MagTense.get_norm_magnetic_flux(H))

        fig = plt.figure()
        ax = fig.gca()
        ax.plot(z_coor, B_h_COMSOL[(id_run[z-1]):id_run[z]], 'bx', label='COMSOL')
        ax.plot(z_coor, MagTense.get_norm_magnetic_flux(H), 'r*', label='MagTense')
        ax.legend()
        ax.set_xlabel('height')
        ax.set_ylabel('B_h')
        plt.title("Iron tiles - " + str(z/100))
        plt.show()

main()

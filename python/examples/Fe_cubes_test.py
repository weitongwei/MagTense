import os
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../source')
import MagTense
import MagTenseStandalone
import util_plot

def plot_H_diff_heights(steps, label):
    fig = plt.figure()
    ax = fig.gca()

    for i in range(10, 11):
        height = i/100
        # Defining grid
        # eval_offset = [0.45, 0.55, height/2]
        eval_offset = [0.65, 0.55, height/2]
        places = [10, 10, 1]
        area = [1, 1, height]
        # Defining occupied places in grid
        filled_positions = [[4, 5, 0],[6, 5, 0]]
        # Defining angle of magnetization in spherical coordinates (azimuth, polar angle) for each tile
        mag_angles = [[math.pi/2, math.pi/2], []] # Magnetization in x-direction with standard value of 1.2 T

        # Optional parameters for setup: n_magnets, filled_positions, mag_angles, eval_points, eval_mode, B_rem
        (tiles, points, grid) = MagTense.setup(places, area, filled_positions=filled_positions, mag_angles=mag_angles, eval_points=[30, 30, 1])
        tiles.set_as_Fe([1])
        tiles.refinement_prism(1)

        struc = np.ones(len(steps))

        if label == 'x':
            pts = np.c_[steps, struc*eval_offset[1], struc*eval_offset[2]]
        elif label == 'y':
            pts = np.c_[struc*eval_offset[0], steps, struc*eval_offset[2]]
        else:
            pts = np.c_[struc*eval_offset[0], struc*eval_offset[1], steps]

        # Standard parameters in settings: max_error=0.00001, max_it=500
        # All steps in one command - iterate_solution = True, return_field = True
        # N is not reused for calculating the H field
        # updated_tiles = MagTense.iterate_magnetization(tiles)
        (updated_tiles, H_field) = MagTense.run_simulation(tiles, pts, grid=grid, plot=False)

        # Plotting B-field for each axis through center
        ax.plot(steps, MagTense.get_norm_magnetic_field(H_field))
        ax.set_xlabel(label)
        #ax.set_ylabel('H_norm [A/m]')
        ax.set_title('H_norm [A/m]')
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

        (updated_tiles, H_field) = MagTense.run_simulation(tiles, points, grid=grid, plot=True)

    plt.grid()
    plt.show()

def main():
    # Set costum points for comparison with COMSOL
    steps = np.linspace(0,1,1000)
    steps_z = np.linspace(-0.5,0.5,1000)

    plot_H_diff_heights(steps, 'x')
    #plot_H_diff_heights(steps, 'y')
    #plot_H_diff_heights(steps_z, 'z')

main()


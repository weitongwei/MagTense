import os
import sys
import numpy as np
import math

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../source')
import MagTense
import MagTenseStandalone
import util_plot

def main():
    # Defining grid
    places = [10, 10, 1]
    area = [1, 1, 0.1]
    # Defining occupied places in grid
    filled_positions = [[3, 3, 0], [3, 5, 0], [5, 3, 0], [5, 5, 0]]
    # Defining angle of magnetization in spherical coordinates (azimuth, polar angle) for each tile
    mag_angles = [[math.pi/2, math.pi/2], [], [], [math.pi/2, 2*math.pi]]

    # Optional parameters for setup: n_magnets, filled_positions, mag_angles, eval_points, eval_mode, B_rem
    (tiles, points, grid) = MagTense.setup(places, area, filled_positions=filled_positions, mag_angles=mag_angles)

    tiles.set_mu_r_ea(1.06)
    tiles.set_mu_r_oa(1.17)

    tiles.set_as_Fe([1,2])

    # Standard parameters in settings: max_error=0.00001, max_it=500
    # All steps in one command - iterate_solution = True, return_field = True
    # N is not reused for calculating the H field
    (updated_tiles, H) = MagTense.run_simulation(tiles, points, grid=grid, plot=True)

main()

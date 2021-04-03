"""Plot a grid of PSFs generated by AIROPA.

Data parameters:
    dir_psf (str) - Folder with the PSF file
    data_psf (str) - PSF file
    data_grid (str) - Grid file

Graphic parameters:
    plots_side (int) - Number of PSF to plots on the side
    buffer (int) - Buffer of the PSF images (px)
"""


from os import path

from astropy.io import fits
from matplotlib import colors, pyplot
import numpy as np


# Data parameters
dir_psf = '/g2/scratch/skterry/work/AIROPA/sean_tests/gc_sky_17/fit/variable'
data_psf = 'c2006_psf_grid.fits'
data_grid = 'c2006_grid_pos.fits'

# Graphic parameters
plots_side = 7
buffer = 130

# Load data
print('Program started')
hdul = fits.open(path.join(dir_psf, data_psf))
psf_data = hdul[0].data
hdul.close()
hdul = fits.open(path.join(dir_psf, data_grid))
grid_data = hdul[0].data
hdul.close()
side_psf = psf_data.shape[1]
n_psfs = psf_data.shape[0]
side_psfs = int(np.sqrt(n_psfs))

# Process data
xs = grid_data[:, 0]
ys = grid_data[:, 1]
min_x = np.min(xs)
max_x = np.max(xs)
min_y = np.min(ys)
max_y = np.max(ys)
grid_x = np.linspace(min_x, max_x, num=plots_side)
grid_y = np.flip(np.linspace(min_y, max_y, num=plots_side))
grid_xs, grid_ys = np.meshgrid(grid_x, grid_y)
out_grid = np.empty((plots_side, plots_side,))
out_grid[:] = np.nan
plot_max = np.max(psf_data)


# Show plots
i_ax = 0
fig, axs = pyplot.subplots(plots_side, plots_side, figsize=(9, 9))

for ax in axs.flat:
    dist = np.hypot((xs - grid_xs.flatten()[i_ax]), (ys - grid_ys.flatten()[i_ax]))
    near = np.argmin(dist)

    ax.imshow(psf_data[near], cmap='jet', aspect='equal', norm=colors.LogNorm(vmin=0.0001, vmax=plot_max))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim([buffer, (side_psf - buffer)])
    ax.set_ylim([buffer, (side_psf - buffer)])
    i_ax += 1


#pyplot.title('TT Star Off-Axis', fontsize=18, x=-3.2, y=8.5)
#pyplot.scatter(1, 1, s=15, color='orange', marker='+')

# End program
pyplot.savefig('test.png', dpi=500)
pyplot.show()
print('Program terminated')

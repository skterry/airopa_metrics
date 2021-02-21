"""Run AIROPA in different modes on a on-sky GC image and produce plots to
compare the astrometry and photometry accuracy.

Data parameters:
    dir_label (string) - 'label.dat' folder
    dir_test (string) - Test data folder
    img_test (str) - Name of the test observation to be plotted

Analysis parameters:
    n_detect_min (int) - Minimum number of detections for a star to be listed
    bright_mag (float) - Magnitude cut to define "bright" stars. Should be where
        residuals in the 'mdr.png' and 'mdm.png' output plots start to increase
        due to noise

Plot parameters:
    img_orig_range (float [2]) - Range of values to show on the original image
    img_res_range (float [2]) - Range of values to show on the residual image
    test_pos (int [8, 2]) - Approximate x and y positions of eight bright test
        stars (2 top-left, 2 top-right, 2 bottom-left, 2 bottom-right of the
        image)
    test_r (float) - Radius of the test stars' imaes (px)
    m_range (float [2]) - Range of magnitudes plotted
    m_range_color (float [2]) - Range of magnitudes colored
    m_histmax (int) - Maximum height in the magnitude histogram
    m_histbin (float) - Bin size in the magnitude histogram
    m_std_range (float [2]) - Range of magnitude STD plotted
    m_std_range_color (float [2]) - Range of magnitude STD colored
    m_std_text (float) - Exponential magnitude STD of the text
    r_std_range (float [2]) - Range of position STD plotted (px)
    r_std_range_color (float [2]) - Range of position STD colored (px)
    r_std_text (float) - Exponential position STD of the text (px)
    fvu_range (float [2]) - Minimum and maximum FVU shown
    fvu_ratio_range (float [2]) - Minimum and maximum FVU ratio shown
    fvu_range_color (float [2]) - Minimum and maximum FVU colored
    fvu_text (float) - Exponential FVU of the text
    quiv_scale (float) - Quiver plot scale
    quiv_legend (float) - Quantity used for the quiver plot legend
    dr_histbin (float) - Bin size in the position difference histogram
    dr_range (float) - Maximum positional difference shown (px)
    dr_histmax (int) - Maximum height in the position difference histogram
    dm_histbin (float) - Bin size in the magnitude difference histogram
    dm_range (float) - Maximum magnitude difference shown (px)
    dm_histmax (int) - Maximum height in the magnitude difference histogram
    m_fake_range (float [2]) - Range of fake stars' magnitudes plotted
    m_fake_range_color (float [2]) - Range of fake stars' magnitudes colored
    m_fake_histmax (int) - Maximum height in the fake stars' magnitude histogram
    m_fake_histbin (float) - Bin size in the fake stars' magnitude histogram
    m_miss_range (float [2]) - Range of missed stars' magnitudes plotted
    m_miss_range_color (float [2]) - Range of missed stars' magnitudes colored
    m_miss_histmax (int) - Maximum height in the missed stars' magnitude
        histogram
    m_miss_histbin (float) - Bin size in the missed stars' magnitude histogram
"""

from airopa_test.analyze_gc import analyze_gc


# Data parameters
dir_label = '/g/lu/data/gc/source_list'
dir_test = ('/g/lu/scratch/jlu/work/ao/airopa/AIROPA_TEST/' +
            'airopa_benchmarks/gc_sky_crash')
img_test = 'c2010'

# Analysis parameters
n_detect_min = 5
bright_mag = 13

# Plot parameters
img_orig_range = [0, 12000]
img_res_range = [-1000, 1000]
test_pos = [[324, 698], [327, 575], [931, 788], [708, 732], [115, 167],
            [68, 287], [599, 111], [931, 428]]
test_r = 50
m_range = [8, 18]
m_range_color = [9, 16]
m_histmax = 200
m_histbin = 0.5
m_std_range = [-2.1, -0.5]
m_std_range_color = [0, 0.1]
m_std_text = -0.8
r_std_range = [-1.4, -0.2]
r_std_range_color = [0, 0.25]
r_std_text = -0.4
fvu_range = [-3.3, 1]
fvu_ratio_range = [0, 2]
fvu_range_color = [-2.5, 0]
fvu_text = 0.5
quiv_scale = 5
quiv_legend = 1
dr_histbin = 0.05
dr_range = 1
dr_histmax = 200
dm_histbin = 0.02
dm_range = 0.3
dm_histmax = 250
m_fake_range = [9, 18]
m_fake_range_color = [12, 17]
m_fake_histbin = 0.5
m_fake_histmax = 60
m_miss_range = [11, 21]
m_miss_range_color = [14, 19]
m_miss_histbin = 0.5
m_miss_histmax = 200

# Start program
analyze_gc(dir_label, dir_test, img_test, n_detect_min, bright_mag,
           img_orig_range, img_res_range, test_pos, test_r, m_range,
           m_range_color, m_std_range, m_std_range_color, m_std_text,
           r_std_range, r_std_range_color, r_std_text, fvu_range,
           fvu_ratio_range, fvu_range_color, quiv_scale, quiv_legend,
           dr_histbin, dr_range, dr_histmax, dm_histbin, dm_range, dm_histmax,
           m_histmax, m_histbin, fvu_text, m_fake_range, m_fake_range_color,
           m_fake_histbin, m_fake_histmax, m_miss_histbin, m_miss_histmax,
           m_miss_range, m_miss_range_color)

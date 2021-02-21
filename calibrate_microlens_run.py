"""Calibrate the AIROPA output catalogs.

Data parameters:
    dir_test (str) - Test data folder
    dir_calib (str) - Photometric calibration folder
"""

from airopa_test.calibrate_microlens import calibrate_microlens


# Data parameters
dir_test = ('/g/lu/scratch/jlu/work/ao/airopa/AIROPA_TEST/airopa_benchmarks/' +
            'gc_sky')
dir_calib = '/g/lu/data/gc/source_list'

# Start program
calibrate_gc(dir_test, dir_calib)

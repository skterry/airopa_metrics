from glob import glob
from os import chdir, path

from nirc2.reduce import calibrate


def calibrate_microlens(dir_test, dir_calib):

    # Start program
    print("\nAnalysis started")

    # Fixed parameters
    modes_str = ['legacy', 'single', 'variable']  # Name of AIROPA modes
    coo_star = 'ob150029'  # Star name in the ".coo" file

    # Prepare the arguments
    args_long = '-f 1 -T 0.0 -I ' + coo_star + ' -N ' + \
                path.join(dir_calib, 'photo_calib.dat') + ' -M K -c 4 '

    for i_mode in modes_str:

        # Find input data
        print()
        print("--- AIROPA \"{0}\" mode ---".format(i_mode))
        cats = glob(path.join(dir_test, 'fit', i_mode, '*_stf.lis'))
        cats = [path.basename(cats[i]) for i in range(len(cats))]
        chdir(path.join(dir_test, 'fit', i_mode))

        # Calibrate
        for i_image in range(len(cats)):
            print("\rCalibrating catalog {0} / {1}...".format((i_image + 1),
                  len(cats)), end="")
            args_temp = args_long + cats[i_image]
            args = args_temp.split()
            calibrate.main(argv=args)

    # End program
    print()
    print()
    print("Analysis completed\n")

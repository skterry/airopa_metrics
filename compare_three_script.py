"""Compare the AIROPA output catalogs.

Parameters:
    old_fit_dir (str) - Output 1 fit/ directory
    new_fit_dir (str) - Output 2 fit/ directory
    old_img_dir (str) - Output 1 image/ directory
    new_img_dir (str) - Output 2 image/ directory
    old_img_root (str) - Output 1 frame
    new_img_root (str) - Output 2 frame
    mode (str) - AIROPA mode
    old_res (str) - Output 1 residual directory
    new_res (str) - Output 2 residual directory
    star_pos (tuple) - Test star position
    star2_pos (tuple) - Test star 2 position (if necessary)
    star3_pos (tuple) - Test star 3 position (if necessary)
"""

from airopa_test.compare_runs import plot_three_psfs,plot_residual_image,compare_fvu_on_three_star,compare_fvu_vs_mag

#Plot two psfs params
one_fit_dir = '../new-POOR-run/fit/'
one_img_dir = '../new-POOR-run/image/'
one_img_root = 'c0271'

two_fit_dir = '../gc_sky_17/fit/'
two_img_dir = '../gc_sky_17/image/'
two_img_root = 'c2006'

three_fit_dir = '../BEST-run/fit/'
three_img_dir = '../BEST-run/image/'
three_img_root = 'c0108'

mode = 'variable'

old_res = one_fit_dir + mode + '/' + one_img_root + '_res.fits'
new_res = two_fit_dir + mode + '/' + two_img_root + '_res.fits'
three_res = three_fit_dir + mode + '/' + three_img_root + '_res.fits'

star_name = 'irs33n'
star_pos = [524,403]
star2_pos = [552,435]
star3_pos = []

#Choose which function to run:
#plot_three_psfs(one_fit_dir, two_fit_dir, three_fit_dir, one_img_root, two_img_root, three_img_root, mode)
#plot_residual_image(old_res,new_res,mode)
compare_fvu_on_three_star(one_fit_dir, two_fit_dir, three_fit_dir, one_img_dir, two_img_dir, three_img_dir, one_img_root, two_img_root, three_img_root, mode, star_name, star_pos, star2_pos, star3_pos)
#compare_fvu_vs_mag(one_fit_dir, two_fit_dir, one_img_dir, two_img_dir, one_img_root, two_img_root, mode)

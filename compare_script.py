"""Calibrate the AIROPA output catalogs.

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
"""

from airopa_test.compare_runs import plot_two_psfs,plot_residual_image,compare_fvu_on_star,compare_fvu_vs_mag

#Plot two psfs params
old_fit_dir = '../new-POOR-run/fit/'
old_img_dir = '../new-POOR-run/image/'
old_img_root = 'c0270'

new_fit_dir = '../new-POOR-run/fit/'
new_img_dir = '../new-POOR-run/image/'
new_img_root = 'c0270'
mode = 'variable'
old_res = old_fit_dir + mode + '/' + old_img_root + '_res.fits'
new_res = new_fit_dir + mode + '/' + new_img_root + '_res.fits'
star_pos = [599, 111]

#Choose which function to run:
#plot_two_psfs(old_fit_dir,new_fit_dir,old_img_root,new_img_root,mode)
#plot_residual_image(old_res,new_res,mode)
compare_fvu_on_star(old_fit_dir, new_fit_dir, old_img_dir, new_img_dir, old_img_root,new_img_root, mode, star_pos)
#compare_fvu_vs_mag(old_fit_dir,new_fit_dir,old_img_dir,new_img_dir,old_img_root,new_img_root,mode)

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from flystar import starlists
from flystar import match
from matplotlib.colors import LogNorm, SymLogNorm
from astropy.nddata import Cutout2D
import pdb

def plot_two_psfs(old_fit_dir, new_fit_dir, old_img_root, new_img_root,mode):
    if mode is 'variable':
        new_psf = new_fit_dir + mode + '/' + new_img_root + '_on_axis_psf.fits'
        old_psf = old_fit_dir + mode + '/' + old_img_root + '_on_axis_psf.fits'
    else:
        new_psf = new_fit_dir + mode + '/' + new_img_root + '_psf.fits'
        old_psf = old_fit_dir + mode + '/' + old_img_root + '_psf.fits'

    npsf = fits.getdata(new_psf)
    opsf = fits.getdata(old_psf)

    min_flux = npsf[npsf > 0].min()
    max_flux = npsf.max() * 0.8
    if min_flux < (1e-5 * max_flux):
        min_flux = 1e-5 * max_flux

    cnorm = LogNorm(vmin=min_flux, vmax=max_flux)
    cmap = 'hot'

    plt.style.use('dark_background')
    plt.figure(figsize=(10, 5))
    plt.subplot(121)
    plt.imshow(opsf, norm=cnorm, cmap=cmap)
    plt.title('Old ' + mode)
    plt.gca().invert_yaxis()

    plt.subplot(122)
    plt.imshow(npsf, norm=cnorm, cmap=cmap)
    plt.title('New ' + mode)
    plt.gca().invert_yaxis()

    plt.show()
    return

def plot_residual_image(old_res, new_res, mode):
    nres = fits.getdata(new_res)
    ores = fits.getdata(old_res)
    
    cmap = 'hot'
    cnorm = SymLogNorm(linthresh=30, linscale=0.5)

    #
    # Zoomed out
    #
    fig = plt.figure(figsize=(10, 5))
    ax1 = plt.subplot(121, aspect='equal')
    im = plt.imshow(ores, cmap=cmap, norm=cnorm)
    plt.title('Old ' + mode)
    plt.gca().invert_yaxis()

    plt.subplot(122, sharex=ax1, sharey=ax1, aspect='equal')
    plt.imshow(nres, cmap=cmap, norm=cnorm)
    plt.title('New ' + mode)
    plt.gca().invert_yaxis()
    
    plt.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)

    #
    # Zoomed in.
    #
    c_pos = [500, 550]
    c_siz = [200, 200]
    nres_c = Cutout2D(nres, c_pos, c_siz)
    ores_c = Cutout2D(ores, c_pos, c_siz)
    cnorm = SymLogNorm(linthresh=30, linscale=0.5, vmin=nres_c.data.min(), vmax=nres_c.data.max())
    
    fig2 = plt.figure(figsize=(10, 5))
    ax2 = plt.subplot(121, aspect='equal')
    im = plt.imshow(ores_c.data, 
                    extent=np.array(ores_c.bbox_original)[::-1].flatten(), 
                    cmap=cmap, norm=cnorm)
    plt.title('Old ' + mode)
    plt.gca().invert_yaxis()

    plt.subplot(122, sharex=ax2, sharey=ax2, aspect='equal')
    plt.imshow(nres_c.data, 
               extent=np.array(nres_c.bbox_original)[::-1].flatten(), 
               cmap=cmap, norm=cnorm)
    plt.title('New ' + mode)
    plt.gca().invert_yaxis()
    
    plt.subplots_adjust(right=0.8)
    cbar_ax = fig2.add_axes([0.85, 0.15, 0.05, 0.7])
    fig2.colorbar(im, cax=cbar_ax)
    
    plt.show()
    return


def compare_fvu_on_star(old_fit_dir, new_fit_dir, old_img_dir, new_img_dir, old_img_root, new_img_root, mode, star_pos=[494, 549]):
    new_res = new_fit_dir + mode + '/' + new_img_root + '_res.fits'
    old_res = old_fit_dir + mode + '/' + old_img_root + '_res.fits'
    img_file = new_img_dir + new_img_root + '.fits'

    nres = fits.getdata(new_res)
    ores = fits.getdata(old_res)
    imag = fits.getdata(img_file)

    c_pos = star_pos
    c_siz = [50, 50]

    # Get the image cutouts
    res_cut_new = Cutout2D(nres, c_pos, c_siz)
    res_cut_old = Cutout2D(ores, c_pos, c_siz)
    img_cut = Cutout2D(imag, c_pos, c_siz)

    # Calculate the FVU for old and new.
    fvu_new = res_cut_new.data.var() / img_cut.data.var()
    fvu_old = res_cut_old.data.var() / img_cut.data.var()

    # Fetch the FVU for this star from the starlists.
    sl_old_name = old_fit_dir + mode + '/' + old_img_root + '_0.8_stf_cal.lis'
    sl_new_name = new_fit_dir + mode + '/' + new_img_root + '_0.8_stf_cal.lis'
    sl_old_mname = old_fit_dir + mode + '/' + old_img_root + '_0.8_metrics.txt'
    sl_new_mname = new_fit_dir + mode + '/' + new_img_root + '_0.8_metrics.txt'    
    sl_old = starlists.StarList.from_lis_file(sl_old_name, fvu_file=sl_old_mname, error=False)
    sl_new = starlists.StarList.from_lis_file(sl_new_name, fvu_file=sl_new_mname, error=False)

    idx_old = np.where((sl_old['x'] > c_pos[0]-5) & (sl_old['x'] < c_pos[0]+5) & 
                       (sl_old['y'] > c_pos[1]-5) & (sl_old['y'] < c_pos[1]+5))[0]
    idx_new = np.where((sl_new['x'] > c_pos[0]-5) & (sl_new['x'] < c_pos[0]+5) & 
                       (sl_new['y'] > c_pos[1]-5) & (sl_new['y'] < c_pos[1]+5))[0]

    print(sl_old[idx_old[0]])
    print(sl_new[idx_new[0]])
    print()

    air_fvu_old = sl_old[idx_old[0]]['fvu']
    air_fvu_new = sl_new[idx_new[0]]['fvu']

    print('    My FVU old = {0:f}, new = {1:f}'.format(fvu_old, fvu_new))
    print('AIROPA FVU old = {0:f}, new = {1:f}'.format(sl_old[idx_old[0]]['fvu'], sl_new[idx_new[0]]['fvu']))

    print('    My FVU ratio: FVU New / Old: {0:f}'.format(fvu_new / fvu_old))
    print('AIROPA FVU ratio: FVU New / Old: {0:f}'.format(sl_new[idx_new[0]]['fvu'] / sl_old[idx_old[0]]['fvu']))
    
    cmap = 'hot'
    vmin_r = np.min([res_cut_old.data.min(), res_cut_new.data.min()])
    vmax_r = np.max([res_cut_old.data.max(), res_cut_new.data.max()])
    vmin_i = img_cut.data.min()
    vmax_i = img_cut.data.max()
    #cnorm_i = SymLogNorm(linthresh=30, linscale=0.5, vmin=img_cut.data.min(), vmax=img_cut.data.max())
    #cnorm_r = SymLogNorm(linthresh=30, linscale=0.5, vmin=vmin_r, vmax=vmax_r)
    
    fig = plt.figure(figsize=(10, 5))
    
    ax = plt.subplot(131, aspect='equal')
    plt.imshow(img_cut.data, 
               extent=np.array(img_cut.bbox_original)[::-1].flatten(), 
               cmap=cmap, vmin=vmin_i, vmax=vmax_i)
    plt.colorbar(orientation='horizontal')
    plt.title('Img ' + mode)

    plt.subplot(132, sharex=ax, sharey=ax, aspect='equal')
    plt.imshow(res_cut_old.data, 
               extent=np.array(res_cut_old.bbox_original)[::-1].flatten(), 
               cmap=cmap, vmin=vmin_r, vmax=vmax_r)
    plt.colorbar(orientation='horizontal')
    plt.setp(plt.gca().get_yticklabels(), visible=False)
    plt.title('Old ' + mode)

    plt.subplot(133, sharex=ax, sharey=ax, aspect='equal')
    plt.imshow(res_cut_new.data, 
               extent=np.array(res_cut_new.bbox_original)[::-1].flatten(), 
               cmap=cmap, vmin=vmin_r, vmax=vmax_r)
    plt.colorbar(orientation='horizontal')
    plt.setp(plt.gca().get_yticklabels(), visible=False)
    plt.title('New ' + mode)
    plt.gca().invert_yaxis()

    plt.show()
    return

def compare_fvu_vs_mag(old_fit_dir, new_fit_dir, old_img_dir, new_img_dir, old_img_root, new_img_root, mode):
    new_res = new_fit_dir + mode + '/' + new_img_root + '_res.fits'
    old_res = old_fit_dir + mode + '/' + old_img_root + '_res.fits'
    img_file = old_img_dir + old_img_root + '.fits'

    nres = fits.getdata(new_res)
    ores = fits.getdata(old_res)
    imag = fits.getdata(img_file)
    #pdb.set_trace()

    # Fetch the FVU for this star from the starlists.
    sl_old_name = old_fit_dir + mode + '/' + old_img_root + '_0.8_stf_cal.lis'
    sl_new_name = new_fit_dir + mode + '/' + new_img_root + '_0.8_stf_cal.lis'
    sl_old_mname = old_fit_dir + mode + '/' + old_img_root + '_0.8_metrics.txt'
    sl_new_mname = new_fit_dir + mode + '/' + new_img_root + '_0.8_metrics.txt'    
    sl_old = starlists.StarList.from_lis_file(sl_old_name, fvu_file=sl_old_mname, error=False)
    sl_new = starlists.StarList.from_lis_file(sl_new_name, fvu_file=sl_new_mname, error=False)

    # Cross match the two starlists.
    idx_old, idx_new, dr, dm = match.match(sl_old['x'], sl_old['y'], sl_old['m'],
                                           sl_new['x'], sl_new['y'], sl_new['m'], dr_tol=30.0, dm_tol=0.5) #Changed dr_tol

    
    air_fvu_old = sl_old[idx_old]['fvu']
    air_fvu_new = sl_new[idx_new]['fvu']
    fvu_old = np.zeros(len(idx_old), dtype=float)
    fvu_new = np.zeros(len(idx_new), dtype=float)
    
    for ii in range(len(idx_old)):
        oo = idx_old[ii]
        nn = idx_new[ii]
        
        c_pos = [sl_new['x'][nn], sl_new['y'][nn]]
        c_siz = [10, 10]

        # Get the image cutouts
        res_cut_new = Cutout2D(nres, c_pos, c_siz)
        res_cut_old = Cutout2D(ores, c_pos, c_siz)
        img_cut = Cutout2D(imag, c_pos, c_siz)

        # Calculate the FVU for old and new.
        fvu_new[ii] = res_cut_new.data.var() / img_cut.data.var()
        fvu_old[ii] = res_cut_old.data.var() / img_cut.data.var()


    plt.figure(figsize=(16,4))
    plt.semilogy(sl_new['m'][idx_new], fvu_old, 'b+', label='Old', alpha=0.5)
    plt.semilogy(sl_new['m'][idx_new], fvu_new, 'rx', label='New', alpha=0.5)
    plt.legend()
    plt.xlabel('Magnitude')
    plt.ylabel('FVU')

    plt.figure(figsize=(16,4))
    plt.plot(sl_new['m'][idx_new], fvu_new / fvu_old, 'k.')
    plt.xlabel('Magnitude')
    plt.ylabel('FVU New / Old')
    plt.ylim(0, 2)
    plt.axhline(1, linestyle='--')
    plt.text(9.5, 0.8, 'New is Better')
    plt.text(9.5, 1.2, 'New is Worse')

    plt.show()
        

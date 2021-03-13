from datetime import datetime
from glob import glob
from os import makedirs, path
import warnings

from astropy import stats, table
from astropy.io import ascii, fits
from astropy.utils.exceptions import AstropyUserWarning
from matplotlib import pyplot, colors, ticker
import pylab as plt
import numpy
import pdb

from flystar import align, starlists



def analyze_microlens(dir_label, dir_test, img_test, n_detect_min, bright_mag,
               img_orig_range, img_res_range, test_pos, test_r, m_range,
               m_range_color, m_std_range, m_std_range_color, m_std_text,
               r_std_range, r_std_range_color, r_std_text, fvu_range,
               fvu_ratio_range, fvu_range_color, quiv_scale, quiv_legend,
               dr_histbin, dr_range, dr_histmax, dm_histbin, dm_range,
               dm_histmax, m_histmax, m_histbin, fvu_text, m_fake_range,
               m_fake_range_color, m_fake_histbin, m_fake_histmax,
               m_miss_histbin, m_miss_histmax, m_miss_range,
               m_miss_range_color):

    # Start program
    print("\nAnalysis started")

    # Fixed parameters
    modes_str = ['legacy', 'single', 'variable']  # Name of AIROPA modes
    xy_buffer = 100  # Buffer in x-y images (px)

    # Load input data
    imgs = glob(path.join(dir_test, 'image', '*.fits'))
    imgs = [path.basename(imgs[i])[0: -5] for i in range(len(imgs))]
    imgs = [x for x in imgs if '_' not in x]
    img_orig = fits.open(path.join(dir_test, 'image', (img_test + '.fits')))
    img_size = [img_orig[0].header['NAXIS1'], img_orig[0].header['NAXIS2']]
    label = table.Table.read(path.join(dir_label, 'label.dat'), format='ascii')
    #pdb.set_trace()
    label.rename_column('col1', 'name')
    label.rename_column('col2', 'm')
    label.rename_column('col3', 'x0')
    label.rename_column('col4', 'y0')
    label.rename_column('col5', 'xerr')
    label.rename_column('col6', 'yerr')
    label.rename_column('col7', 'vx')
    label.rename_column('col8', 'vy')
    label.rename_column('col9', 'vxerr')
    label.rename_column('col10', 'vyerr')
    label.rename_column('col11', 't')
    label.remove_column('col12')
    label.remove_column('col13')
    label['x0'] /= -0.009952
    label['y0'] /= 0.009952
    label['xerr'] /= 0.009952
    label['yerr'] /= 0.009952
    label['vx'] /= -9.952
    label['vy'] /= 9.952
    label['vxerr'] /= 9.952
    label['vyerr'] /= 9.952
    idx_label = numpy.where(label['name'] == 'ob150029')[0][0]
    coo = table.Table.read(path.join(dir_test, 'image', (img_test + '.coo')),
                           format='ascii')
    dx_label = coo['col1'][0] - label['x0'][idx_label]
    dy_label = coo['col2'][0] - label['y0'][idx_label]
    label['x0'] += dx_label
    label['y0'] += dy_label
    cat = starlists.StarList.from_lis_file(path.join(dir_test, 'fit', 'legacy',
                                           (img_test + '_0.8_stf_cal.lis')),
                                           error=False)
    label['dt'] = cat['t'][0] - label['t']
    label['x'] = label['x0'] + (label['vx'] * label['dt'])
    label['y'] = label['y0'] + (label['vy'] * label['dt'])
    sl_label_all = starlists.StarList.from_table(label)
    sl_label = sl_label_all[numpy.all([(sl_label_all['x'] > 0),
                                       (sl_label_all['x'] < 1024),
                                       (sl_label_all['y'] > 0),
                                       (sl_label_all['y'] < 1024)], axis=0)]

    # Save the parameters
    folder_out = path.join(dir_test, 'plot')

    if not path.exists(folder_out):
        makedirs(folder_out)

    now = datetime.now()
    params = open(path.join(dir_test, 'plot', 'params_plt.txt'), 'w')
    params.write('Date analyzed:            ' + str(now.year) + '-' +
                 str(now.month) + '-' + str(now.day) + '\n\n')
    params.write('--- PARAMETERS ---' + '\n')
    params.write('Number of images:         ' + str(len(imgs)) + '\n')
    params.write('Test image:               ' + img_test + '\n')
    params.write('Minimum detection number: ' + str(n_detect_min) + '\n')
    params.close()

    # Find test stars
    cat_test = glob(path.join(path.join(dir_test, 'fit', 'variable'),
                              img_test + '_*_stf.lis'))
    tab_test = table.Table.read(cat_test[0], format='ascii', delimiter='\s')
    tab_test.rename_column('col4', 'x')
    tab_test.rename_column('col5', 'y')
    tab_test.rename_column('col2', 'm')

    for i_test in range(8):
        test_dist = numpy.hypot((tab_test['x'] - test_pos[i_test][0]),
                                (tab_test['y'] - test_pos[i_test][1]))
        test_idx = numpy.argmin(test_dist)
        test_pos[i_test] = [int(tab_test['x'][test_idx]),
                            int(tab_test['y'][test_idx])]

    # Show original image
    img_ps = 0.009942  # Pixel scale
    img_ticks = 2
    img_shape = img_orig[0].data.shape
    xy_px = numpy.array(img_shape)
    xy_ticks_max = numpy.floor((xy_px * img_ps) / (2 * img_ticks)) * img_ticks
    x_ticks = numpy.arange(-xy_ticks_max[1], (xy_ticks_max[1] + img_ticks), img_ticks)
    y_ticks = numpy.arange(-xy_ticks_max[0], (xy_ticks_max[0] + img_ticks), img_ticks)
    x_ticks_orig = (x_ticks / img_ps) + ((xy_px[1] - 1) / 2)
    y_ticks_orig = (y_ticks / img_ps) + ((xy_px[0] - 1) / 2)
    x_ticks_str = [str(i) for i in numpy.around(x_ticks, decimals=10)]
    y_ticks_str = [str(i) for i in numpy.around(y_ticks, decimals=10)]
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    plt.imshow(img_orig[0].data, clim=img_orig_range, cmap='hot',
                  aspect='equal')
    ax.set_position([0.18, 0.1, 0.75, 0.75])
    ax.invert_yaxis()
    ax.set_xticks(x_ticks_orig)
    ax.set_yticks(y_ticks_orig)
    ax.set_xticklabels(x_ticks_str)
    ax.set_yticklabels(y_ticks_str)
    plt.xlabel("x (\")", fontsize=18)
    plt.ylabel("y (\")", fontsize=18)
    plt.tick_params(labelsize=15)
    save_path = path.join(folder_out, 'image')
    plt.savefig(save_path + '.png')
    plt.savefig(save_path + '.pdf')
    plt.close(fig)

    # Show original test stars
    test_str = ['tl1', 'tl2', 'tr1', 'tr2', 'bl1', 'bl2', 'br1', 'br2']
    ms = [None] * 3
    fvus = [None] * 3
    img_ticks_zoom = 0.2
    xy_ticks_max_zoom = numpy.floor((xy_px * img_ps) / (2 * img_ticks_zoom)) * img_ticks
    x_ticks_zoom = numpy.arange(-xy_ticks_max_zoom[1], (xy_ticks_max_zoom[1] + img_ticks_zoom), img_ticks_zoom)
    y_ticks_zoom = numpy.arange(-xy_ticks_max_zoom[0], (xy_ticks_max_zoom[0] + img_ticks_zoom), img_ticks_zoom)
    x_ticks_orig_zoom = (x_ticks_zoom / img_ps) + ((xy_px[1] - 1) / 2)
    y_ticks_orig_zoom = (y_ticks_zoom / img_ps) + ((xy_px[0] - 1) / 2)
    x_ticks_str_zoom = [str(i) for i in numpy.around(x_ticks_zoom, decimals=10)]
    y_ticks_str_zoom = [str(i) for i in numpy.around(y_ticks_zoom, decimals=10)]

    for i_test in range(8):
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)
        plt.imshow(img_orig[0].data, clim=img_orig_range, cmap='hot',
                      aspect='equal')
        ax.set_xticks(x_ticks_orig_zoom)
        ax.set_yticks(y_ticks_orig_zoom)
        ax.set_xticklabels(x_ticks_str_zoom)
        ax.set_yticklabels(y_ticks_str_zoom)
        plt.axis([(test_pos[i_test][0] - test_r),
                     (test_pos[i_test][0] + test_r),
                     (test_pos[i_test][1] - test_r),
                     (test_pos[i_test][1] + test_r)])
        ax.set_position([0.2, 0.1, 0.7, 0.7])
        plt.xlabel("x (\")", fontsize=18)
        plt.ylabel("y (\")", fontsize=18)
        plt.tick_params(labelsize=15)
        save_path = path.join(folder_out, ('image_test_' + test_str[i_test]))
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

    # Analyze AIROPA modes
    sts = [None] * 3
    sts_detect = [None] * 3
    trs = [None] * 3
    stsls = []
    sl_calib = []

    for i_mode in range(0, 3):

        # Load output data
        print()
        print("--- AIROPA \"{0}\" mode ---".format(modes_str[i_mode]))
        folder_in = path.join(dir_test, 'fit', modes_str[i_mode])
        img_res = fits.open(path.join(folder_in, (img_test + '_res.fits')))
        cat_out = starlists.StarList.from_lis_file(
            glob(path.join(folder_in, (img_test + '_*_stf.lis')))[0],
            error=False)
        metrics_out_test = \
            ((table.Table.read(glob(path.join(folder_in,(img_test +
                                                         '_*_metrics.txt')))[0],
                               format='ascii', data_start=0)).columns[0]) ** 2

        # Analyze test observation
        bright_idx_test = numpy.where(cat_out['m'] <= bright_mag)[0]
        fvu_bright_test = metrics_out_test[bright_idx_test]
        fvu_bright_filt_test = \
            stats.sigma_clip(numpy.ma.masked_invalid(fvu_bright_test), sigma=3,
                             maxiters=5)
        fvu_bright_mean_test = numpy.mean(fvu_bright_filt_test)

        # Analyze observations
        sls = []
        sl_test = None

        for i_image, name_image in enumerate(imgs):

            # Load output data
            print("\rLoading image {0} / {1}...".format((i_image + 1),
                  len(imgs)), end="")
            cat_out1 = glob(path.join(folder_in, (name_image + '_*_stf.lis')))
            sls.append(starlists.StarList.from_lis_file(cat_out1[0],
                       error=False))
            cat_out2 = glob(path.join(folder_in, (name_image +
                            '_*_metrics.txt')))
            metrics_out = table.Table.read(cat_out2[0], format='ascii',
                                           data_start=0)
            sls[-1].add_column(table.Column(metrics_out[metrics_out.colnames[0]]**2), name='fvu')
            sls[-1].add_column(table.Column(numpy.hypot((sls[-1]['x'] - (img_size[0] / 2)),
                                                        (sls[-1]['y'] - (img_size[1] / 2)))), name='r')

            if name_image == img_test:
                cat_out3 = glob(path.join(folder_in, (name_image +
                                '_*_stf_cal.lis')))
                calib = table.Table.read(path.join(folder_in, cat_out3[0]),
                                         format='ascii')
                irs_stars = numpy.where([i[0].isdigit() for i in
                                        calib['col1']])[0]

                for i_irs in irs_stars:
                    if len(calib['col1'][i_irs]) < 6:
                        calib['col1'][i_irs] = 'irs' + calib['col1'][i_irs]

                calib.rename_column('col1', 'name')
                calib.rename_column('col2', 'm')
                calib.rename_column('col3', 't')
                calib.rename_column('col4', 'x')
                calib.rename_column('col5', 'y')
                calib.rename_column('col6', 'snr')
                calib.rename_column('col7', 'corr')
                calib.rename_column('col8', 'N_frames')
                calib.rename_column('col9', 'flux')
                sl_calib.append(starlists.StarList.from_table(calib))
                sl_test = sls[-1]

        print()

        # Match stars
        print("Matching stars...")
        warnings.simplefilter('ignore')
        sts[i_mode], trs[i_mode] = align.mosaic_lists(sls,
                                                      trans_args=[{'order': 0},
                                                                  {'order': 0}],
                                                      update_ref_per_iter=False,
                                                      verbose=False)
        sts[i_mode].combine_lists('fvu')
        sts[i_mode].combine_lists('r')
        detect_idx = numpy.where(sts[i_mode]['n_detect'] > n_detect_min)[0]
        sts_detect[i_mode] = sts[i_mode][detect_idx]
        name = sts_detect[i_mode]['name']
        n_detect = sts_detect[i_mode]['n_detect']
        m = sts_detect[i_mode]['m_avg']
        ms[i_mode] = m
        m_std = sts_detect[i_mode]['m_std']
        fvu = sts_detect[i_mode]['fvu_avg']
        fvus[i_mode] = fvu
        r = sts_detect[i_mode]['r_avg']
        x = sts_detect[i_mode]['x_avg']
        n_stars = len(x)
        y = sts_detect[i_mode]['y_avg']
        r_std = (sts_detect[i_mode]['x_std'] + sts_detect[i_mode]['y_std']) / 2
        bright_idx = numpy.where(sts_detect[i_mode]['m_avg'] <= bright_mag)[0]
        m_bright_std = m_std[bright_idx]
        r_bright_std = r_std[bright_idx]
        fvu_bright = fvu[bright_idx]
        # warnings.filterwarnings('ignore', category=RuntimeWarning)
        m_bright_std_filt =\
            stats.sigma_clip(numpy.ma.masked_invalid(m_bright_std), sigma=3,
                             maxiters=5)
        r_bright_std_filt =\
            stats.sigma_clip(numpy.ma.masked_invalid(r_bright_std), sigma=3,
                             maxiters=5)
        fvu_bright_filt = \
            stats.sigma_clip(numpy.ma.masked_invalid(fvu_bright), sigma=3,
                             maxiters=5)
        # warnings.filterwarnings('default', category=RuntimeWarning)
        m_bright_std_mean = numpy.mean(m_bright_std_filt)
        r_bright_std_mean = numpy.mean(r_bright_std_filt)
        fvu_bright_mean = numpy.mean(fvu_bright_filt)
        print()
        print(numpy.median(fvu_bright_filt))
        print()
        stsl = starlists.StarList.from_table(table.Table([name, x, y, r_std, m,
                                             m_std, r, fvu, n_detect],
                                             names=('name', 'x', 'y',
                                                    'xy_std_orig', 'm',
                                                    'm_std_orig', 'r', 'fvu',
                                                    'n_detect_orig')))
        stsls.append(stsl)
        calib_label = [sl_calib[i_mode], sl_label]
        calib_match, _ = align.mosaic_lists(calib_label,
                                            dr_tol=[2, 2],
                                            dm_tol=[0.5, 0.5],
                                            trans_args=[{'order': 3},
                                                        {'order': 3}],
                                            init_guess_mode='name',
                                            update_ref_per_iter=False,
                                            verbose=False)
        idx_fake = numpy.where([numpy.any(calib_match['detect'][i] == 0) and not
                               numpy.any(calib_match['detect'][i] == 1) for i in
                               range(len(calib_match))])[0]
        idx_miss = numpy.where([numpy.any(calib_match['detect'][i] == 1) and not
                               numpy.any(calib_match['detect'][i] == 0) for i in
                               range(len(calib_match))])[0]
        name_fake = [calib_match['name_in_list'][i][0] for i in idx_fake]
        name_miss = [calib_match['name_in_list'][i][1] for i in idx_miss]
        idx_fake2 = [numpy.where(sl_calib[i_mode]['name'] == i)[0][0] for i in
                     name_fake]
        idx_miss2 = [numpy.where(sl_label['name'] == i)[0][0] for i in
                     name_miss]
        sl_fake = sl_test[idx_fake2]
        sl_miss = sl_label[idx_miss2]

        # Show residual image
        print("Producing plots...")
        folder_out = path.join(dir_test, 'plot', modes_str[i_mode])

        if not path.exists(folder_out):
            makedirs(folder_out)

        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)
        plt.imshow(img_res[0].data, clim=img_res_range, cmap='hot',
                      aspect='equal')
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        ax.invert_yaxis()
        plt.xlabel("x (\")", fontsize=18)
        plt.ylabel("y (\")", fontsize=18)
        plt.tick_params(labelsize=15)
        save_path = path.join(folder_out, 'image_res')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Show residual test stars
        for i_test in range(8):
            fig = plt.figure(figsize=(6, 6))
            ax = fig.add_subplot(111)
            plt.imshow(img_res[0].data, clim=img_res_range, cmap='hot',
                          aspect='equal')
            ax.set_xticks(x_ticks_orig_zoom)
            ax.set_yticks(y_ticks_orig_zoom)
            ax.set_xticklabels(x_ticks_str_zoom)
            ax.set_yticklabels(y_ticks_str_zoom)
            plt.scatter((cat_out['x'] - 1), (cat_out['y'] - 1), marker='o',
                           c='g', s=70)
            plt.axis([(test_pos[i_test][0] - test_r),
                         (test_pos[i_test][0] + test_r),
                         (test_pos[i_test][1] - test_r),
                         (test_pos[i_test][1] + test_r)])
            ax.set_position([0.2, 0.1, 0.7, 0.7])
            plt.xlabel("x (\")", fontsize=18)
            plt.ylabel("y (\")", fontsize=18)
            plt.tick_params(labelsize=15)
            save_path = path.join(folder_out, ('image_res_test_' +
                                  test_str[i_test]))
            plt.savefig(save_path + '.png')
            plt.savefig(save_path + '.pdf')
            plt.close(fig)

        # Magnitude vs magnitude STD
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        plt.scatter(m, m_std, s=10, c='k')
        ax.set_yscale('log')
        plt.xlim(m_range)
        plt.ylim([(10 ** i) for i in m_std_range])
        plt.xlabel("m", fontsize=18)
        plt.ylabel(r"m$_{\sigma}$", fontsize=18)
        plt.tick_params(axis='y', which='minor')
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f"))
        ax.yaxis.set_minor_formatter(ticker.FormatStrFormatter("%.2f"))
        plt.tick_params(which='both', labelsize=12)
        save_path = path.join(folder_out, 'mmstd')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Magnitude vs magnitude STD with text
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        plt.scatter(m, m_std, s=10, c='k')
        plt.text((((m_range[1] - m_range[0]) / 15) + m_range[0]),
                 (10 ** m_std_text),
                 (r"m$_{\sigma,m\leq" + str(bright_mag) + "}$="
                  + " {0:.4f}".format(m_bright_std_mean)), ha='left',
                 fontsize=15)
        ax.set_yscale('log')
        plt.xlim(m_range)
        plt.ylim([(10 ** i) for i in m_std_range])
        plt.xlabel("m", fontsize=18)
        plt.ylabel(r"m$_{\sigma}$", fontsize=18)
        plt.tick_params(labelsize=15)
        save_path = path.join(folder_out, 'mmstd_text')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Magnitude vs magnitude STD, colored by position
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        scatter = plt.scatter(m, m_std, s=10, c=r,
                                 cmap=plt.get_cmap('plasma'), vmin=0,
                                 vmax=(numpy.hypot(img_size[0], img_size[1]) /
                                       2))
        ax.set_yscale('log')
        plt.xlim(m_range)
        plt.ylim([(10 ** i) for i in m_std_range])
        plt.xlabel("m", fontsize=18)
        plt.ylabel(r"m$_{\sigma}$", fontsize=18)
        plt.tick_params(labelsize=15)
        colorbar_ax = fig.add_axes([0.18, 0.9, 0.75, 0.03])
        fig.colorbar(scatter, cax=colorbar_ax, orientation='horizontal')
        colorbar_ax.xaxis.set_ticks_position('top')
        colorbar_ax.tick_params(labelsize=15)
        colorbar_ax.set_xlabel("R (px)", fontsize=18)
        save_path = path.join(folder_out, 'mmstdr')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Position vs magnitude STD
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        plt.scatter(r, m_std, s=10, c='k')
        ax.set_yscale('log')
        plt.xlim([0, (numpy.hypot(img_size[0], img_size[1]) / 2)])
        plt.ylim([(10 ** i) for i in m_std_range])
        plt.xlabel("R (px)", fontsize=18)
        plt.ylabel(r"m$_{\sigma}$", fontsize=18)
        plt.tick_params(labelsize=15)
        save_path = path.join(folder_out, 'rmstd')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Position vs magnitude STD, colored by magnitude
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        scatter = plt.scatter(r, m_std, s=10, c=m, cmap=plt.get_cmap('plasma'),
                              vmin=m_range_color[0], vmax=m_range_color[1])
        ax.set_yscale('log')
        plt.xlim([0, (numpy.hypot(img_size[0], img_size[1]) / 2)])
        plt.ylim([(10 ** i) for i in m_std_range])
        plt.xlabel("R (px)", fontsize=18)
        plt.ylabel(r"m$_{\sigma}$", fontsize=18)
        plt.tick_params(labelsize=15)
        colorbar_ax = fig.add_axes([0.18, 0.9, 0.75, 0.03])
        fig.colorbar(scatter, cax=colorbar_ax, orientation='horizontal')
        colorbar_ax.xaxis.set_ticks_position('top')
        colorbar_ax.tick_params(labelsize=15)
        colorbar_ax.set_xlabel("m", fontsize=18)
        save_path = path.join(folder_out, 'rmstdm')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Magnitude STD map for bright stars
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        scatter = plt.scatter(x, y, s=10, c=m_std,
                                 cmap=plt.get_cmap('plasma'),
                                 vmin=m_std_range_color[0],
                                 vmax=m_std_range_color[1])
        plt.xlim([-xy_buffer, (img_size[0] + xy_buffer)])
        plt.ylim([-xy_buffer, (img_size[1] + xy_buffer)])
        ax.set_aspect('equal', 'box')
        plt.xlabel("x (px)", fontsize=18)
        plt.ylabel("y (px)", fontsize=18)
        plt.tick_params(labelsize=15)
        colorbar_ax = fig.add_axes([0.18, 0.9, 0.75, 0.03])
        fig.colorbar(scatter, cax=colorbar_ax, orientation='horizontal')
        colorbar_ax.xaxis.set_ticks_position('top')
        colorbar_ax.tick_params(labelsize=15)
        colorbar_ax.set_xlabel(r"m$_{\sigma}$", fontsize=18)
        save_path = path.join(folder_out, 'xymstd')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Magnitude vs position STD
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        plt.scatter(m, r_std * img_ps * 1e3, s=10, c='k')
        ax.set_yscale('log')
        plt.xlim(m_range)
        plt.ylim([(10 ** i) * img_ps * 1e3 for i in r_std_range])
        plt.xlabel("m", fontsize=18)
        plt.ylabel(r"r$_{\sigma}$ (mas)", fontsize=18)
        plt.tick_params(axis='y', which='minor')
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
        ax.yaxis.set_minor_formatter(ticker.FormatStrFormatter("%.1f"))
        plt.tick_params(which='both', labelsize=12)
        save_path = path.join(folder_out, 'mrstd')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Magnitude vs position STD with text
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        plt.scatter(m, r_std, s=10, c='k')
        plt.text((((m_range[1] - m_range[0]) / 15) + m_range[0]),(10 **
                    r_std_text), (r"r$_{\sigma,m\leq" + str(bright_mag) + "}$="
                    + " {0:.4f} mas".format(r_bright_std_mean * img_ps * 1e3)), ha='left',
                    fontsize=15)
        ax.set_yscale('log')
        plt.xlim(m_range)
        plt.ylim([(10 ** i) for i in r_std_range])
        plt.xlabel("m", fontsize=18)
        plt.ylabel(r"r$_{\sigma}$ (px)", fontsize=18)
        plt.tick_params(labelsize=15)
        save_path = path.join(folder_out, 'mrstd_text')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Magnitude vs position STD, colored by position
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        scatter = plt.scatter(m, r_std, s=10, c=r,
                                 cmap=plt.get_cmap('plasma'), vmin=0,
                                 vmax=(numpy.hypot(img_size[0], img_size[1]) /
                                       2))
        ax.set_yscale('log')
        plt.xlim(m_range)
        plt.ylim([(10 ** i) for i in r_std_range])
        plt.xlabel("m", fontsize=18)
        plt.ylabel(r"r$_{\sigma}$ (px)", fontsize=18)
        plt.tick_params(labelsize=15)
        colorbar_ax = fig.add_axes([0.18, 0.9, 0.75, 0.03])
        fig.colorbar(scatter, cax=colorbar_ax, orientation='horizontal')
        colorbar_ax.xaxis.set_ticks_position('top')
        colorbar_ax.tick_params(labelsize=15)
        colorbar_ax.set_xlabel("R (px)", fontsize=18)
        save_path = path.join(folder_out, 'mrstdr')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Position vs position STD
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        plt.scatter(r, r_std, s=10, c='k')
        ax.set_yscale('log')
        plt.xlim([0, (numpy.hypot(img_size[0], img_size[1]) / 2)])
        plt.ylim([(10 ** i) for i in r_std_range])
        plt.xlabel("R (px)", fontsize=18)
        plt.ylabel(r"r$_{\sigma}$ (px)", fontsize=18)
        plt.tick_params(labelsize=15)
        save_path = path.join(folder_out, 'rrstd')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Position vs position STD, colored by magnitude
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        scatter = plt.scatter(r, r_std, s=10, c=m,
                                 cmap=plt.get_cmap('plasma'),
                                 vmin=m_range_color[0], vmax=m_range_color[1])
        ax.set_yscale('log')
        plt.xlim([0, (numpy.hypot(img_size[0], img_size[1]) / 2)])
        plt.ylim([(10 ** i) for i in r_std_range])
        plt.xlabel("R (px)", fontsize=18)
        plt.ylabel(r"r$_{\sigma}$ (px)", fontsize=18)
        plt.tick_params(labelsize=15)
        colorbar_ax = fig.add_axes([0.18, 0.9, 0.75, 0.03])
        fig.colorbar(scatter, cax=colorbar_ax, orientation='horizontal')
        colorbar_ax.xaxis.set_ticks_position('top')
        colorbar_ax.tick_params(labelsize=15)
        colorbar_ax.set_xlabel("m", fontsize=18)
        save_path = path.join(folder_out, 'rrstdm')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Position STD map for bright stars
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        scatter = plt.scatter(x, y, s=10, c=r_std,
                                 cmap=plt.get_cmap('plasma'),
                                 vmin=r_std_range_color[0],
                                 vmax=r_std_range_color[1])
        plt.xlim([-xy_buffer, (img_size[0] + xy_buffer)])
        plt.ylim([-xy_buffer, (img_size[1] + xy_buffer)])
        ax.set_aspect('equal', 'box')
        plt.xlabel("x (px)", fontsize=18)
        plt.ylabel("y (px)", fontsize=18)
        plt.tick_params(labelsize=15)
        colorbar_ax = fig.add_axes([0.18, 0.9, 0.75, 0.03])
        fig.colorbar(scatter, cax=colorbar_ax, orientation='horizontal')
        colorbar_ax.xaxis.set_ticks_position('top')
        colorbar_ax.tick_params(labelsize=15)
        colorbar_ax.set_xlabel(r"r$_{\sigma}$ (px)", fontsize=18)
        save_path = path.join(folder_out, 'xyrstd')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Magnitude vs FVU
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        plt.scatter(m, fvu, s=10, c='k')
        ax.set_yscale('log')
        plt.xlim(m_range)
        plt.ylim([(10 ** i) for i in fvu_range])
        plt.xlabel("m", fontsize=18)
        plt.ylabel("FVU", fontsize=18)
        plt.tick_params(labelsize=15)
        save_path = path.join(folder_out, 'mfvu')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Magnitude vs FVU with text
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        plt.scatter(m, fvu, s=10, c='k')
        plt.text((((m_range[1] - m_range[0]) / 15) + m_range[0]), (10 **
                                                                      fvu_text),
                    (r"m$_{FVU,m\leq" + str(bright_mag) + "}$=" +
                     " {0:.4f}".format(fvu_bright_mean)), ha='left',
                    fontsize=15)
        ax.set_yscale('log')
        plt.xlim(m_range)
        plt.ylim([(10 ** i) for i in fvu_range])
        plt.xlabel("m", fontsize=18)
        plt.ylabel("FVU", fontsize=18)
        plt.tick_params(labelsize=15)
        save_path = path.join(folder_out, 'mfvu_text')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Magnitude vs FVU, colored by position
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        scatter = plt.scatter(m, fvu, s=10, c=r,
                                 cmap=plt.get_cmap('plasma'), vmin=0,
                                 vmax=(numpy.hypot(img_size[0], img_size[1]) /
                                       2))
        ax.set_yscale('log')
        plt.xlim(m_range)
        plt.ylim([(10 ** i) for i in fvu_range])
        plt.xlabel("m", fontsize=18)
        plt.ylabel("FVU", fontsize=18)
        plt.tick_params(labelsize=15)
        colorbar_ax = fig.add_axes([0.18, 0.9, 0.75, 0.03])
        fig.colorbar(scatter, cax=colorbar_ax, orientation='horizontal')
        colorbar_ax.xaxis.set_ticks_position('top')
        colorbar_ax.tick_params(labelsize=15)
        colorbar_ax.set_xlabel("R (px)", fontsize=18)
        save_path = path.join(folder_out, 'mfvur')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Position vs FVU
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        plt.scatter(r, fvu, s=10, c='k')
        ax.set_yscale('log')
        plt.xlim([0, (numpy.hypot(img_size[0], img_size[1]) / 2)])
        plt.ylim([(10 ** i) for i in fvu_range])
        plt.xlabel("R (px)", fontsize=18)
        plt.ylabel("FVU", fontsize=18)
        plt.tick_params(labelsize=15)
        save_path = path.join(folder_out, 'rfvu')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Position vs FVU, colored by magnitude
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        scatter = plt.scatter(r, fvu, s=10, c=m,
                                 cmap=plt.get_cmap('plasma'),
                                 vmin=m_range_color[0], vmax=m_range_color[1])
        ax.set_yscale('log')
        plt.xlim([0, (numpy.hypot(img_size[0], img_size[1]) / 2)])
        plt.ylim([(10 ** i) for i in fvu_range])
        plt.xlabel("R (px)", fontsize=18)
        plt.ylabel("FVU", fontsize=18)
        plt.tick_params(labelsize=15)
        colorbar_ax = fig.add_axes([0.18, 0.9, 0.75, 0.03])
        fig.colorbar(scatter, cax=colorbar_ax, orientation='horizontal')
        colorbar_ax.xaxis.set_ticks_position('top')
        colorbar_ax.tick_params(labelsize=15)
        colorbar_ax.set_xlabel("m", fontsize=18)
        save_path = path.join(folder_out, 'rfvum')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # FVU map for bright stars
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        scatter = plt.scatter(x, y, s=10, c=fvu,
                                 cmap=plt.get_cmap('plasma_r'),
                                 norm=colors.LogNorm(
                                     vmin=(10 ** fvu_range_color[0]),
                                     vmax=(10 ** fvu_range_color[1])))
        plt.xlim([-xy_buffer, (img_size[0] + xy_buffer)])
        plt.ylim([-xy_buffer, (img_size[1] + xy_buffer)])
        ax.set_aspect('equal', 'box')
        plt.xlabel("x (px)", fontsize=18)
        plt.ylabel("y (px)", fontsize=18)
        plt.tick_params(labelsize=15)
        colorbar_ax = fig.add_axes([0.18, 0.9, 0.75, 0.03])
        fig.colorbar(scatter, cax=colorbar_ax, orientation='horizontal')
        colorbar_ax.xaxis.set_ticks_position('top')
        colorbar_ax.tick_params(labelsize=15)
        colorbar_ax.set_xlabel(r"FVU", fontsize=18)
        save_path = path.join(folder_out, 'xyfvu')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Position STD vs FVU
        fig = pyplot.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        pyplot.scatter(r_std, fvu, s=10)
        ax.set_xscale('log')
        ax.set_yscale('log')
        pyplot.xlim([(10 ** i) for i in r_std_range])
        pyplot.ylim([(10 ** i) for i in fvu_range])
        pyplot.xlabel(r"r$_{\sigma}$ (px)", fontsize=18)
        pyplot.ylabel("FVU", fontsize=18)
        pyplot.tick_params(labelsize=15)
        save_path = path.join(folder_out, 'rstdfvu')
        pyplot.savefig(save_path + '.png')
        pyplot.savefig(save_path + '.pdf')
        pyplot.close(fig)

        # Magnitude STD vs FVU
        fig = pyplot.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        pyplot.scatter(m_std, fvu, s=10)
        ax.set_xscale('log')
        ax.set_yscale('log')
        pyplot.xlim([(10 ** i) for i in m_std_range])
        pyplot.ylim([(10 ** i) for i in fvu_range])
        pyplot.xlabel(r"m$_{\sigma}$ (px)", fontsize=18)
        pyplot.ylabel("FVU", fontsize=18)
        pyplot.tick_params(labelsize=15)
        save_path = path.join(folder_out, 'mstdfvu')
        pyplot.savefig(save_path + '.png')
        pyplot.savefig(save_path + '.pdf')
        pyplot.close(fig)

        # Histogram of detections by magnitude
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.15, 0.75, 0.75])
        plt.hist(m, bins=numpy.arange(m_range[0], m_range[1], m_histbin),
                    color='b')
        plt.xlim(m_range)
        plt.ylim([0, m_histmax])
        plt.xlabel("m", fontsize=18)
        plt.ylabel("N", fontsize=18)
        plt.tick_params(labelsize=15)
        save_path = path.join(folder_out, 'm')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Histogram of detections by magnitude with text
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.15, 0.75, 0.75])
        plt.hist(m, bins=numpy.arange(m_range[0], m_range[1], m_histbin),
                    color='b')
        plt.text((m_range[0] + 0.3), (3 * m_histmax / 4),
                    (r"$N_{match}$=" + " {0:d}".format(n_stars)), ha='left',
                    fontsize=15)
        plt.xlim(m_range)
        plt.ylim([0, m_histmax])
        plt.xlabel("m", fontsize=18)
        plt.ylabel("N", fontsize=18)
        plt.tick_params(labelsize=15)
        save_path = path.join(folder_out, 'm_text')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Fake stars positions, colored by magnitude
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        scatter = plt.scatter(sl_fake['x'], sl_fake['y'], s=10,
                                 c=sl_fake['m'], cmap=plt.get_cmap('plasma'),
                                 vmin=m_fake_range_color[0],
                                 vmax=m_fake_range_color[1])
        plt.xlim([-xy_buffer, (img_size[0] + xy_buffer)])
        plt.ylim([-xy_buffer, (img_size[1] + xy_buffer)])
        plt.xlabel("x (px)", fontsize=18)
        plt.ylabel("y (px)", fontsize=18)
        plt.tick_params(labelsize=15)
        colorbar_ax = fig.add_axes([0.18, 0.9, 0.75, 0.03])
        fig.colorbar(scatter, cax=colorbar_ax, orientation='horizontal')
        colorbar_ax.xaxis.set_ticks_position('top')
        colorbar_ax.tick_params(labelsize=15)
        colorbar_ax.set_xlabel("m", fontsize=18)
        save_path = path.join(folder_out, 'xym_fake')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Fake stars histogram
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.15, 0.75, 0.75])
        plt.hist(sl_fake['m'], bins=numpy.arange(m_fake_range[0],
                    m_fake_range[1], m_fake_histbin), color='b')
        plt.xlim(m_fake_range)
        plt.ylim([0, m_fake_histmax])
        plt.xlabel("m", fontsize=18)
        plt.ylabel("N", fontsize=18)
        plt.tick_params(labelsize=15)
        save_path = path.join(folder_out, 'm_fake')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Fake stars histogram with text
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.15, 0.75, 0.75])
        plt.hist(sl_fake['m'], bins=numpy.arange(m_fake_range[0],
                    m_fake_range[1], m_fake_histbin), color='b')
        plt.text((m_fake_range[0] + 0.3), (3 * m_fake_histmax / 4),
                    (r"$N_{fake}$=" + " {0:d}".format(len(sl_fake))), ha='left',
                    fontsize=15)
        plt.xlim(m_fake_range)
        plt.ylim([0, m_fake_histmax])
        plt.xlabel("m", fontsize=18)
        plt.ylabel("N", fontsize=18)
        plt.tick_params(labelsize=15)
        save_path = path.join(folder_out, 'm_fake_text')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Missed stars positions, colored by magnitude
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        scatter = plt.scatter(sl_miss['x'], sl_miss['y'], s=10,
                                 c=sl_miss['m'], cmap=plt.get_cmap('plasma'),
                                 vmin=m_miss_range_color[0],
                                 vmax=m_miss_range_color[1])
        plt.xlim([-xy_buffer, (img_size[0] + xy_buffer)])
        plt.ylim([-xy_buffer, (img_size[1] + xy_buffer)])
        plt.xlabel("x (px)", fontsize=18)
        plt.ylabel("y (px)", fontsize=18)
        plt.tick_params(labelsize=15)
        colorbar_ax = fig.add_axes([0.18, 0.9, 0.75, 0.03])
        fig.colorbar(scatter, cax=colorbar_ax, orientation='horizontal')
        colorbar_ax.xaxis.set_ticks_position('top')
        colorbar_ax.tick_params(labelsize=15)
        colorbar_ax.set_xlabel("m", fontsize=18)
        save_path = path.join(folder_out, 'xym_miss')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Missed stars histogram
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.15, 0.75, 0.75])
        plt.hist(sl_miss['m'], bins=numpy.arange(m_miss_range[0],
                    m_miss_range[1], m_miss_histbin), color='b')
        plt.xlim(m_miss_range)
        plt.ylim([0, m_miss_histmax])
        plt.xlabel("m", fontsize=18)
        plt.ylabel("N", fontsize=18)
        plt.tick_params(labelsize=15)
        save_path = path.join(folder_out, 'm_miss')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Missed stars histogram with text
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.15, 0.75, 0.75])
        plt.hist(sl_miss['m'], bins=numpy.arange(m_miss_range[0],
                    m_miss_range[1], m_miss_histbin), color='b')
        plt.text((m_miss_range[0] + 0.3), (3 * m_miss_histmax / 4),
                    (r"$N_{miss}$=" + " {0:d}".format(len(sl_miss))), ha='left',
                    fontsize=15)
        plt.xlim(m_miss_range)
        plt.ylim([0, m_miss_histmax])
        plt.xlabel("m", fontsize=18)
        plt.ylabel("N", fontsize=18)
        plt.tick_params(labelsize=15)
        save_path = path.join(folder_out, 'm_miss_text')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Test image magnitude vs FVU with text
        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(111)
        ax.set_position([0.18, 0.1, 0.75, 0.75])
        plt.scatter(cat_out['m'], metrics_out_test, s=10, c='k')
        plt.text((((m_range[1] - m_range[0]) / 15) + m_range[0]), (10 **
                                                                   fvu_text),
                 (r"m$_{FVU,m\leq" + str(bright_mag) + "}$=" +
                  " {0:.4f}".format(fvu_bright_mean_test)), ha='left',
                 fontsize=15)
        ax.set_yscale('log')
        plt.xlim(m_range)
        plt.ylim([(10 ** i) for i in fvu_range])
        plt.xlabel("m", fontsize=18)
        plt.ylabel("FVU", fontsize=18)
        plt.tick_params(labelsize=15)
        save_path = path.join(folder_out, 'mfvu_test_text')
        plt.savefig(save_path + '.png')
        plt.savefig(save_path + '.pdf')
        plt.close(fig)

        # Save catalog
        st_txt = sts_detect[i_mode]
        st_txt.combine_lists('flux')
        st_txt.remove_rows(numpy.where(st_txt['flux_avg'] == -1))
        ascii.write(st_txt, output=path.join(folder_out, 'cat.txt'),
                    include_names=['x_avg', 'y_avg', 'flux_avg'],
                    overwrite=True)
        lines = open(path.join(folder_out, 'cat.txt')).readlines()
        with open(path.join(folder_out, 'cat.txt'), 'w') as file:
            file.writelines(lines[1:])

    # Match modes
    print()
    print("Matching modes...")
    # warnings.filterwarnings('ignore', category=AstropyUserWarning)
    
    st, tr = align.mosaic_lists(stsls, trans_args=[{'order': 0}, {'order': 0}],
                                update_ref_per_iter=False, verbose=False)
    # warnings.filterwarnings('default', category=AstropyUserWarning)
    st_all = st[numpy.where(st['n_detect'] == 3)]
    st_all.add_column((st_all['x'][:, 1] - st_all['x'][:, 0]), name='dx_01')
    st_all.add_column((st_all['y'][:, 1] - st_all['y'][:, 0]), name='dy_01')
    st_all.add_column(numpy.hypot((st_all['x'][:, 1] - st_all['x'][:, 0]),
                      (st_all['y'][:, 1] - st_all['y'][:, 0])), name='dr_01')
    st_all.add_column((st_all['m'][:, 1] - st_all['m'][:, 0]), name='dm_01')
    st_all.add_column((st_all['fvu'][:, 1] / st_all['fvu'][:, 0]),
                      name='fvur_01')
    st_all.add_column((st_all['x'][:, 2] - st_all['x'][:, 0]), name='dx_02')
    st_all.add_column((st_all['y'][:, 2] - st_all['y'][:, 0]), name='dy_02')
    st_all.add_column(numpy.hypot((st_all['x'][:, 2] - st_all['x'][:, 0]),
                      (st_all['y'][:, 2] - st_all['y'][:, 0])), name='dr_02')
    st_all.add_column((st_all['m'][:, 2] - st_all['m'][:, 0]), name='dm_02')
    st_all.add_column((st_all['fvu'][:, 2] / st_all['fvu'][:, 0]),
                      name='fvur_02')
    st_all.add_column((st_all['x'][:, 2] - st_all['x'][:, 1]), name='dx_12')
    st_all.add_column((st_all['y'][:, 2] - st_all['y'][:, 1]), name='dy_12')
    st_all.add_column(numpy.hypot((st_all['x'][:, 2] - st_all['x'][:, 1]),
                      (st_all['y'][:, 2] - st_all['y'][:, 1])), name='dr_12')
    st_all.add_column((st_all['m'][:, 2] - st_all['m'][:, 1]), name='dm_12')
    st_all.add_column((st_all['fvu'][:, 2] / st_all['fvu'][:, 1]),
                      name='fvur_12')

    # Quiver plot of position difference, legacy vs single
    print("Producing plots...")
    folder_out = path.join(dir_test, 'plot')
    fig = plt.figure(figsize=(6, 6.5))
    ax = fig.add_subplot(111)
    ax.set_position([0.18, 0.1, 0.75, 0.75])
    quiv = plt.quiver(st_all['x'][:, 0], st_all['y'][:, 0], st_all['dx_01'],
                         st_all['dy_01'], scale=quiv_scale,
                         scale_units='height', width=0.003)
    plt.quiverkey(quiv, 0.4, 0.9, quiv_legend, (str(quiv_legend) + " px"),
                     labelpos='W', coordinates='figure', fontproperties={'size':
                     18})
    plt.xlim([-xy_buffer, (img_size[0] + xy_buffer)])
    plt.ylim([-xy_buffer, (img_size[1] + xy_buffer)])
    ax.set_aspect('equal', 'box')
    plt.xlabel("x (px)", fontsize=18)
    plt.ylabel("y (px)", fontsize=18)
    plt.tick_params(labelsize=15)
    save_path = path.join(folder_out, 'xydxdy_01')
    plt.savefig(save_path + '.png')
    plt.savefig(save_path + '.pdf')
    plt.close(fig)

    # Quiver plot of position difference, legacy vs variable
    fig = plt.figure(figsize=(6, 6.5))
    ax = fig.add_subplot(111)
    ax.set_position([0.18, 0.1, 0.75, 0.75])
    quiv = plt.quiver(st_all['x'][:, 0], st_all['y'][:, 0], st_all['dx_02'],
                         st_all['dy_02'], scale=quiv_scale,
                         scale_units='height', width=0.003)
    plt.quiverkey(quiv, 0.4, 0.9, quiv_legend, (str(quiv_legend) + " px"),
                     labelpos='W', coordinates='figure', fontproperties={'size':
                     18})
    plt.xlim([-xy_buffer, (img_size[0] + xy_buffer)])
    plt.ylim([-xy_buffer, (img_size[1] + xy_buffer)])
    ax.set_aspect('equal', 'box')
    plt.xlabel("x (px)", fontsize=18)
    plt.ylabel("y (px)", fontsize=18)
    plt.tick_params(labelsize=15)
    save_path = path.join(folder_out, 'xydxdy_02')
    plt.savefig(save_path + '.png')
    plt.savefig(save_path + '.pdf')
    plt.close(fig)

    # Quiver plot of position difference, single vs variable
    fig = plt.figure(figsize=(6, 6.5))
    ax = fig.add_subplot(111)
    ax.set_position([0.18, 0.1, 0.75, 0.75])
    quiv = plt.quiver(st_all['x'][:, 1], st_all['y'][:, 1], st_all['dx_12'],
                         st_all['dy_12'], scale=quiv_scale,
                         scale_units='height', width=0.003)
    plt.quiverkey(quiv, 0.4, 0.9, quiv_legend, (str(quiv_legend) + " px"),
                     labelpos='W', coordinates='figure', fontproperties={'size':
                     18})
    plt.xlim([-xy_buffer, (img_size[0] + xy_buffer)])
    plt.ylim([-xy_buffer, (img_size[1] + xy_buffer)])
    ax.set_aspect('equal', 'box')
    plt.xlabel("x (px)", fontsize=18)
    plt.ylabel("y (px)", fontsize=18)
    plt.tick_params(labelsize=15)
    save_path = path.join(folder_out, 'xydxdy_12')
    plt.savefig(save_path + '.png')
    plt.savefig(save_path + '.pdf')
    plt.close(fig)

    # Histogram of position difference, legacy vs single
    fig = plt.figure(figsize=(6, 6.5))
    ax = fig.add_subplot(111)
    ax.set_position([0.18, 0.15, 0.75, 0.75])
    plt.hist(st_all['dr_01'], bins=numpy.arange(0, max(st_all['dr_01']),
                dr_histbin), color='b')
    plt.xlim([0, dr_range])
    plt.ylim([0, dr_histmax])
    plt.xlabel(r"$\Delta$r (px)", fontsize=18)
    plt.ylabel("N", fontsize=18)
    plt.tick_params(labelsize=15)
    save_path = path.join(folder_out, 'dr_01')
    plt.savefig(save_path + '.png')
    plt.savefig(save_path + '.pdf')
    plt.close(fig)

    # Histogram of position difference, legacy vs variable
    fig = plt.figure(figsize=(6, 6.5))
    ax = fig.add_subplot(111)
    ax.set_position([0.18, 0.15, 0.75, 0.75])
    plt.hist(st_all['dr_02'], bins=numpy.arange(0, max(st_all['dr_02']),
                dr_histbin), color='b')
    plt.xlim([0, dr_range])
    plt.ylim([0, dr_histmax])
    plt.xlabel(r"$\Delta$r (px)", fontsize=18)
    plt.ylabel("N", fontsize=18)
    plt.tick_params(labelsize=15)
    save_path = path.join(folder_out, 'dr_02')
    plt.savefig(save_path + '.png')
    plt.savefig(save_path + '.pdf')
    plt.close(fig)

    # Histogram of position difference, single vs variable
    fig = plt.figure(figsize=(6, 6.5))
    ax = fig.add_subplot(111)
    ax.set_position([0.18, 0.15, 0.75, 0.75])
    plt.hist(st_all['dr_12'], bins=numpy.arange(0, max(st_all['dr_12']),
                dr_histbin), color='b')
    plt.xlim([0, dr_range])
    plt.ylim([0, dr_histmax])
    plt.xlabel(r"$\Delta$r (px)", fontsize=18)
    plt.ylabel("N", fontsize=18)
    plt.tick_params(labelsize=15)
    save_path = path.join(folder_out, 'dr_12')
    plt.savefig(save_path + '.png')
    plt.savefig(save_path + '.pdf')
    plt.close(fig)

    # Plot of magnitude difference, legacy vs single
    fig = plt.figure(figsize=(6, 6.5))
    ax = fig.add_subplot(111)
    ax.set_position([0.18, 0.1, 0.75, 0.75])
    scatter = plt.scatter(st_all['x'][:, 0], st_all['y'][:, 0], s=10,
                             c=st_all['dm_01'],
                             cmap=plt.get_cmap('plasma_r'), vmin=-dm_range,
                             vmax=dm_range)
    plt.xlim([-xy_buffer, (img_size[0] + xy_buffer)])
    plt.ylim([-xy_buffer, (img_size[1] + xy_buffer)])
    ax.set_aspect('equal', 'box')
    plt.xlabel("x (px)", fontsize=18)
    plt.ylabel("y (px)", fontsize=18)
    plt.tick_params(labelsize=15)
    colorbar_ax = fig.add_axes([0.18, 0.9, 0.75, 0.03])
    fig.colorbar(scatter, cax=colorbar_ax, orientation='horizontal')
    colorbar_ax.xaxis.set_ticks_position('top')
    colorbar_ax.tick_params(labelsize=15)
    colorbar_ax.set_xlabel(r"$\Delta$m", fontsize=18)
    save_path = path.join(folder_out, 'xydm_01')
    plt.savefig(save_path + '.png')
    plt.savefig(save_path + '.pdf')
    plt.close(fig)

    # Plot of magnitude difference, legacy vs variable
    fig = plt.figure(figsize=(6, 6.5))
    ax = fig.add_subplot(111)
    ax.set_position([0.18, 0.1, 0.75, 0.75])
    scatter = plt.scatter(st_all['x'][:, 0], st_all['y'][:, 0], s=10,
                             c=st_all['dm_02'],
                             cmap=plt.get_cmap('plasma_r'), vmin=-dm_range,
                             vmax=dm_range)
    plt.xlim([-xy_buffer, (img_size[0] + xy_buffer)])
    plt.ylim([-xy_buffer, (img_size[1] + xy_buffer)])
    ax.set_aspect('equal', 'box')
    plt.xlabel("x (px)", fontsize=18)
    plt.ylabel("y (px)", fontsize=18)
    plt.tick_params(labelsize=15)
    colorbar_ax = fig.add_axes([0.18, 0.9, 0.75, 0.03])
    fig.colorbar(scatter, cax=colorbar_ax, orientation='horizontal')
    colorbar_ax.xaxis.set_ticks_position('top')
    colorbar_ax.tick_params(labelsize=15)
    colorbar_ax.set_xlabel(r"$\Delta$m", fontsize=18)
    save_path = path.join(folder_out, 'xydm_02')
    plt.savefig(save_path + '.png')
    plt.savefig(save_path + '.pdf')
    plt.close(fig)

    # Plot of magnitude difference, single vs variable
    fig = plt.figure(figsize=(6, 6.5))
    ax = fig.add_subplot(111)
    ax.set_position([0.18, 0.1, 0.75, 0.75])
    scatter = plt.scatter(st_all['x'][:, 1], st_all['y'][:, 1], s=10,
                             c=st_all['dm_12'],
                             cmap=plt.get_cmap('plasma_r'), vmin=-dm_range,
                             vmax=dm_range)
    plt.xlim([-xy_buffer, (img_size[0] + xy_buffer)])
    plt.ylim([-xy_buffer, (img_size[1] + xy_buffer)])
    ax.set_aspect('equal', 'box')
    plt.xlabel("x (px)", fontsize=18)
    plt.ylabel("y (px)", fontsize=18)
    plt.tick_params(labelsize=15)
    colorbar_ax = fig.add_axes([0.18, 0.9, 0.75, 0.03])
    fig.colorbar(scatter, cax=colorbar_ax, orientation='horizontal')
    colorbar_ax.xaxis.set_ticks_position('top')
    colorbar_ax.tick_params(labelsize=15)
    colorbar_ax.set_xlabel(r"$\Delta$m", fontsize=18)
    save_path = path.join(folder_out, 'xydm_12')
    plt.savefig(save_path + '.png')
    plt.savefig(save_path + '.pdf')
    plt.close(fig)

    # Histogram of magnitude difference, legacy vs single
    fig = plt.figure(figsize=(6, 6.5))
    ax = fig.add_subplot(111)
    ax.set_position([0.18, 0.15, 0.75, 0.75])
    plt.hist(numpy.abs(st_all['dm_01']), bins=numpy.arange(0,
                max(numpy.abs(st_all['dm_01'])), dm_histbin), color='b')
    plt.xlim([0, dm_range])
    plt.ylim([0, dm_histmax])
    plt.xlabel(r"|$\Delta$m|", fontsize=18)
    plt.ylabel("N", fontsize=18)
    plt.tick_params(labelsize=15)
    save_path = path.join(folder_out, 'dm_01')
    plt.savefig(save_path + '.png')
    plt.savefig(save_path + '.pdf')
    plt.close(fig)

    # Histogram of magnitude difference, legacy vs variable
    fig = plt.figure(figsize=(6, 6.5))
    ax = fig.add_subplot(111)
    ax.set_position([0.18, 0.15, 0.75, 0.75])
    plt.hist(numpy.abs(st_all['dm_02']), bins=numpy.arange(0,
                max(numpy.abs(st_all['dm_02'])), dm_histbin), color='b')
    plt.xlim([0, dm_range])
    plt.ylim([0, dm_histmax])
    plt.xlabel(r"|$\Delta$m|", fontsize=18)
    plt.ylabel("N", fontsize=18)
    plt.tick_params(labelsize=15)
    save_path = path.join(folder_out, 'dm_02')
    plt.savefig(save_path + '.png')
    plt.savefig(save_path + '.pdf')
    plt.close(fig)

    # Histogram of magnitude difference, single vs variable
    fig = plt.figure(figsize=(6, 6.5))
    ax = fig.add_subplot(111)
    ax.set_position([0.18, 0.15, 0.75, 0.75])
    plt.hist(numpy.abs(st_all['dm_12']), bins=numpy.arange(0,
                max(numpy.abs(st_all['dm_12'])), dm_histbin), color='b')
    plt.xlim([0, dm_range])
    plt.ylim([0, dm_histmax])
    plt.xlabel(r"|$\Delta$m|", fontsize=18)
    plt.ylabel("N", fontsize=18)
    plt.tick_params(labelsize=15)
    save_path = path.join(folder_out, 'dm_12')
    plt.savefig(save_path + '.png')
    plt.savefig(save_path + '.pdf')
    plt.close(fig)

    # Magnitude vs FVU, legacy vs single
    fig = plt.figure(figsize=(6, 6.5))
    ax = fig.add_subplot(111)
    ax.set_position([0.18, 0.1, 0.75, 0.75])
    plt.scatter(ms[0], fvus[0], s=10, c='b', label=modes_str[0])
    plt.scatter(ms[1], fvus[1], s=10, c='r', label=modes_str[1])
    ax.set_yscale('log')
    plt.xlim(m_range)
    plt.ylim([(10 ** i) for i in fvu_range])
    ax.legend(loc='upper left', fontsize=15)
    plt.xlabel("m", fontsize=18)
    plt.ylabel("FVU", fontsize=18)
    plt.tick_params(labelsize=15)
    save_path = path.join(folder_out, 'mfvu_01')
    plt.savefig(save_path + '.png')
    plt.savefig(save_path + '.pdf')
    plt.close(fig)

    # Magnitude vs FVU, legacy vs variable
    fig = plt.figure(figsize=(6, 6.5))
    ax = fig.add_subplot(111)
    ax.set_position([0.18, 0.1, 0.75, 0.75])
    plt.scatter(ms[0], fvus[0], s=10, c='b', label=modes_str[0])
    plt.scatter(ms[2], fvus[2], s=10, c='r', label=modes_str[2])
    ax.set_yscale('log')
    plt.xlim(m_range)
    plt.ylim([(10 ** i) for i in fvu_range])
    ax.legend(loc='upper left', fontsize=15)
    plt.xlabel("m", fontsize=18)
    plt.ylabel("FVU", fontsize=18)
    plt.tick_params(labelsize=15)
    save_path = path.join(folder_out, 'mfvu_02')
    plt.savefig(save_path + '.png')
    plt.savefig(save_path + '.pdf')
    plt.close(fig)

    # Magnitude vs FVU, single vs variable
    fig = plt.figure(figsize=(6, 6.5))
    ax = fig.add_subplot(111)
    ax.set_position([0.18, 0.1, 0.75, 0.75])
    plt.scatter(ms[1], fvus[1], s=10, c='b', label=modes_str[1])
    plt.scatter(ms[2], fvus[2], s=10, c='r', label=modes_str[2])
    ax.set_yscale('log')
    plt.xlim(m_range)
    plt.ylim([(10 ** i) for i in fvu_range])
    ax.legend(loc='upper left', fontsize=15)
    plt.xlabel("m", fontsize=18)
    plt.ylabel("FVU", fontsize=18)
    plt.tick_params(labelsize=15)
    save_path = path.join(folder_out, 'mfvu_12')
    plt.savefig(save_path + '.png')
    plt.savefig(save_path + '.pdf')
    plt.close(fig)

    # Magnitude vs FVU difference, legacy vs single
    fig = plt.figure(figsize=(6, 6.5))
    ax = fig.add_subplot(111)
    ax.set_position([0.18, 0.1, 0.75, 0.75])
    plt.scatter(st_all['m'][:, 0], st_all['fvur_01'], s=10, c='k')
    plt.plot(m_range, [1, 1], 'r--')
    plt.xlim(m_range)
    plt.ylim(fvu_ratio_range)
    plt.xlabel("m", fontsize=18)
    plt.ylabel(r"FVU$_{single}$ / FVU$_{legacy}$", fontsize=18)
    plt.tick_params(labelsize=15)
    save_path = path.join(folder_out, 'mfvur_01')
    plt.savefig(save_path + '.png')
    plt.savefig(save_path + '.pdf')
    plt.close(fig)

    # Magnitude vs FVU difference, legacy vs variable
    fig = plt.figure(figsize=(6, 6.5))
    ax = fig.add_subplot(111)
    ax.set_position([0.18, 0.1, 0.75, 0.75])
    plt.scatter(st_all['m'][:, 0], st_all['fvur_02'], s=10, c='k')
    plt.plot(m_range, [1, 1], 'r--')
    plt.xlim(m_range)
    plt.ylim(fvu_ratio_range)
    plt.xlabel("m", fontsize=18)
    plt.ylabel(r"FVU$_{variable}$ / FVU$_{legacy}$", fontsize=18)
    plt.tick_params(labelsize=15)
    save_path = path.join(folder_out, 'mfvur_02')
    plt.savefig(save_path + '.png')
    plt.savefig(save_path + '.pdf')
    plt.close(fig)

    # Magnitude vs FVU difference, single vs variable
    fig = plt.figure(figsize=(6, 6.5))
    ax = fig.add_subplot(111)
    ax.set_position([0.18, 0.1, 0.75, 0.75])
    plt.scatter(st_all['m'][:, 0], st_all['fvur_12'], s=10, c='k')
    plt.plot(m_range, [1, 1], 'r--')
    plt.xlim(m_range)
    plt.ylim(fvu_ratio_range)
    plt.xlabel("m", fontsize=18)
    plt.ylabel(r"FVU$_{variable}$ / FVU$_{single}$", fontsize=18)
    plt.tick_params(labelsize=15)
    save_path = path.join(folder_out, 'mfvur_12')
    plt.savefig(save_path + '.png')
    plt.savefig(save_path + '.pdf')
    plt.close(fig)

    # End program
    print()
    print("Analysis completed\n")
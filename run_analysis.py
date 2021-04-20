import os
import shutil
import numpy as np
from astropy import table
from datetime import date
from astropy.io import ascii
from nirc2.reduce import prep_analysis
from nirc2.reduce import util
from microlens.jlu import analysis, align_flystar
from astroquery.mast import Observations, Catalogs 
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
from astroquery.vsa import Vsa
from astropy import units as u
from flystar import starlists, align, match
from flystar import transforms
#from flystar import analysis, plots
from importlib import reload
import pylab as plt

#------------------------------------------------------------
# User Input
#------------------------------------------------------------

src_list_dir = "/g2/scratch/skterry/work/mb19284/source_list/"
combo_file = "mag20jun25os_mb19284_kp_tdOpen.fits"
epoch = '20jun25os' #YearMonthDay ('os' == osiris)
instrument = 'osiris' #'nirc2' or 'osiris' only two options.
filt = 'kp' #nirc2 or osiris filters (i.e. kp, ks, kn3, j, js, h)
corr = 0.9 #corrMain correleation value

target = 'mb19284'
ra = '18:05:55.06'
dec = '-30:20:12.9'
coo = (1126.90, 1017.74)

#------------------------------------------------------------
# Begin Functions
#------------------------------------------------------------
#def psflist():
#    psfStars = prep_analysis.generate_list("../combo/" + combo_file)
#    x = []
#    finalTable = table.Table()
#    for i in zip(psfStars["x"], psfStars["y"], psfStars["m"]):
#        x.append((list(i)))
#    print(x)
#    return x


#psfStars = [
#   [1007.42, 1068.14, 1],
#   [1205.73, 678.51, 1],
#   [467.15, 1013.42, 1],
#   [1103.31, 615.17, 1],
#   [1865.72, 920.42, 1],
#   [534.72, 1211.68, 1],
#   [1180.87, 1689.58, 1],
#   [1156.67, 583.0, 0]
#            ]

if combo_file.endswith('.fits'):
    short_fname = combo_file[13:-5] #for proper filename reading.

# Generate PSF list and prep Starfinder.
def prepstf():
    finalTable = prep_analysis.generate_list("../combo/" + combo_file)
    psfStars = []
    for i in zip(finalTable["x"], finalTable["y"], finalTable["m"]):
        psfStars.append((list(i)))
    
    prep_analysis.prepStarfinder(src_list_dir, target, psfStars[0], psfStars, filt, instrument)
    
    util.mkdir('../clean/' + short_fname)
    shutil.copyfile('../clean/c.lis', '../clean/' + short_fname + '/c.lis')
    shutil.copyfile('../clean/strehl_source.txt', '../clean/' + short_fname + '/strehl_source.txt')

# Run Starfinder
def runstf():
    analysisObject = analysis.mb19284(epoch, filt)
    analysisObject.starfinderCombo()

if combo_file.endswith('.fits'):
    combo_file_stripped = combo_file[:-5] #for proper filename reading.

# Prep photometric calibration
def prepcalib():
    target_coords = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
    search_rad = 0.25 * u.arcmin

    sourceTable = Vsa.query_region(target_coords, radius=search_rad, programme_id='VVV', database='VVVDR4')
    #print(sourceTable.colnames)
    #print(sourceTable[0:5])
    #print(sourceTable['ra', 'dec', 'ksAperMag1', 'ksAperMag1Err'][0:5])

    targTable = Vsa.query_region(target_coords, radius= 1 * u.arcsec, programme_id='VVV', database='VVVDR4')
    x = targTable["ra"][0]
    y = targTable["dec"][0]
    v4_coords = SkyCoord(x, y, unit=(u.deg, u.deg), frame='icrs')
    #print(v4_coords)
    #print(target_coords)
    d_dec = sourceTable['dec']*u.degree - v4_coords.dec
    #print(d_dec[0:5])
    #print(d_dec.arcsec[0:5])
    d_ras = (sourceTable['ra']*u.degree - v4_coords.ra) * np.cos(v4_coords.dec.radian)
    #print(d_ras[0:5])
    #print(d_ras.arcsec[0:5])
    sourceTable['d_dec'] = d_dec.to_value('arcsec')
    sourceTable['d_ras'] = d_ras.to_value('arcsec')

    # Plot the NIRC2 image with the label.dat list on top of it.
    nirc2_label_file = '../combo/starfinder/' + combo_file_stripped + '_' + str(corr) + '_stf.lis'
    nirc2_label = table.Table.read(nirc2_label_file, format = "ascii", names = ("name",
                                                                            "m",
                                                                            "t0",
                                                                            "Xpix",
                                                                            "Ypix",
                                                                            "snr?",
                                                                            "corr?",
                                                                            "frames?",
                                                                            "flux?"))
    print(nirc2_label)
    nirc2_label["x"] = np.round((np.array(nirc2_label["Xpix"][:]) - (nirc2_label["Xpix"][0])) * (-9.942 / 1000), 4)
    nirc2_label["y"] = np.round((np.array(nirc2_label["Ypix"][:]) - (nirc2_label["Ypix"][0])) * (9.942 / 1000), 4)
    plotter_list = nirc2_label["name", "m", "x", "y"]
    print("***********************************************")
    print(plotter_list)

    prep_analysis.plot_starlist_on_image_arcsec(plotter_list, 
                                          '../combo/' + combo_file, 
                                          coo,
                                          flip = True, 
                                          label = False, 
                                          magCut = 14, 
                                          verbose = False)
    plt.xlim(-10, 10)
    plt.ylim(-10, 10)
    plt.plot(sourceTable["d_ras"], sourceTable["d_dec"], "kx")

    gdx = np.where(nirc2_label['m'] < 13)[0]
    nirc2_label_g = nirc2_label[gdx]
    ndx_g, ndx_n, dr, dm = match.match(sourceTable["d_ras"], sourceTable["d_dec"], sourceTable["ksAperMag3"],
                                   nirc2_label_g["x"], nirc2_label_g["y"], nirc2_label_g["m"],
                                  dr_tol = 0.1, dm_tol = 10)

    nirc2_match = nirc2_label_g[ndx_n]
    sourceTable_match = sourceTable[ndx_g]

    #print(ndx_g)
    #print(nirc2_match['name', 'x', 'y', 'm'])
    def replaceStarWithS00(lst):
        lnew = []
        for s in lst:
            if s.find("star") >= 0:
                s = s.replace("star_", "")
                s_new = "S"
                for i in range(3 - len(s)):
                    s_new += "0"
                lnew.append(s_new + s)
            else:
                 lnew.append(s)
        return lnew

    nirc2_match["name"] = replaceStarWithS00(nirc2_match["name"][:])
    nirc2_match["name"]
    nirc2_calib_file = '/g/lu/data/microlens/source_list/' + "ob120169" + '_photo.dat'
    calib = table.Table.read(nirc2_calib_file, format='ascii.commented_header', delimiter='\s', header_start=-1)
    #print(calib)

    print('****************************')
    print('* TEMPLATE PHOTO.DAT TABLE *')
    print('****************************')

    # Flag to track whether VVVdr4 is in this file.
    found_vvv_dr4 = False

    for cc in calib.meta['comments']:
        print(cc)
    calib.pprint(max_width=-1)

    today = date.today()
    dt = today.strftime("%Y_%m_%d")
    calib.meta["comments"][0] = "# Version: " + dt
    calib.meta['comments'][1] = '# Analysis: point to this exact .ipynb analysis' # Replace this!

    # Getting the reference stars (i.e. non-target matched stars) and 
    # joining them in a string to replace the Filt -- ... part
    not_targ = np.where(nirc2_match['name'] != target)[0]
    ref_stars_str = ','.join(nirc2_match['name'][not_targ])

    # In this case, the filter is Ks_VVV (since we are using VVV to calibrate). So create this column.
    calib['Ks_VVV']= 0.0
    calib['Ks_VVV_err'] = 0.0

    # Create new comment in header with refstars for this filter
    new_comm = 'Filt -- Ks_VVV: Ks band (VVVdr4) -- ' + ref_stars_str

    # Remove other filter comments and add this one
    for i in range(4):
        calib.meta['comments'].pop()
    calib.meta['comments'].append(new_comm)

    # Remove template stars from table
    calib.remove_rows([0, 1, 2, 3])

    # Add the new stars, and set their magnitudes and errors to the Ks_VVV data
    for ss in range(len(nirc2_match)):
        sdx = np.where(calib['Star'] == nirc2_match['name'][ss])[0]
    
        if len(sdx) == 1:
            star = sdx[0]
        
            calib['Ks_VVV'][star] = sourceTable_match['ksAperMag1'][ss]
            calib['Ks_VVV_err'][star] = sourceTable_match['ksAperMag1Err'][ss]
        else:
            calib.add_row()
            star = -1
            calib['Star'][star] = nirc2_match['name'][ss]
            calib['x_pos'][star] = nirc2_match['x'][ss]
            calib['y_pos'][star] = nirc2_match['y'][ss]
            calib['Ks_VVV'][star] = sourceTable_match['ksAperMag1'][ss]
            calib['Ks_VVV_err'][star] = sourceTable_match['ksAperMag1Err'][ss]
        
        
    calib['Star'].format = '{:<13s}'
    calib['x_pos'].format = '{:6.3f}'
    calib['y_pos'].format = '{:6.3f}'
    calib['x_vel'].format = '{:6.3f}'
    calib['y_vel'].format = '{:6.3f}'
    calib['t0'].format = '{:7.2f}'
    calib['var?'].format = '{:1d}'

    calib['J'].format = '{:6.3f}'
    calib['H'].format = '{:6.3f}'
    calib['K'].format = '{:6.3f}'
    calib['Ks_VVV'].format = '{:6.3f}'

    calib['J_err'].format = '{:5.3f}'
    calib['H_err'].format = '{:5.3f}'
    calib['K_err'].format = '{:5.3f}'
    calib['Ks_VVV_err'].format = '{:5.3f}'
    
    print('*************')
    print('* NEW TABLE *')
    print('*************')

    for cc in calib.meta['comments']:
        print(cc)
    calib.pprint(max_width=-1)

    calib_file = src_list_dir + target + '_photo.dat'


    _out = open(calib_file, 'w+')
    for cc in calib.meta['comments']:
        if cc.startswith("#"):
            _out.write('#' + cc + '\n')
        else:
            _out.write('# ' + cc + '\n')


    del calib.meta['comments']
    _out.write('#')
    ascii.write(calib, _out, format='fixed_width', delimiter=' ', bookend=False)
    _out.close()


# Run photometric calibration, generates *_photo.dat file.
def calibrate():
    analysisObject = analysis.mb19284(epoch, filt)
    analysisObject.calibrateCombo()


# Run label, generates *_label.dat file.
def label():
    calStfLisPath = '../combo/starfinder/' + combo_file_stripped + '_' + str(corr) + '_stf_cal.lis'
    calStfLis = align_flystar.starlists.read_starlist(calStfLisPath, error=False)

    calStfLis["xarc"] = (calStfLis["x"] - calStfLis["x"][0]) * (9.942 / 1000) # This is relative distance to target in x
    calStfLis["yarc"] = (calStfLis["y"] - calStfLis["y"][0]) * (9.942 / 1000)  # This is relative distance to target in y
    calStfLis["r2d"] = np.round(np.sqrt((calStfLis["xarc"] ** 2) +  (calStfLis["yarc"] ** 2)), 4) # Euclidean distance
    calStfLis["xarc"] = np.round(calStfLis["xarc"], 4)
    calStfLis["yarc"] = np.round(calStfLis["yarc"], 4)

    calStfLis["xerr"] = 0.0
    calStfLis["yerr"] = 0.0
    calStfLis["vx"] = 0.0
    calStfLis["vy"] = 0.0
    calStfLis["vxerr"] = 0.0
    calStfLis["vyerr"] = 0.0
    calStfLis["use"] = 1
    calStfLis.rename_column("m", "kp")
    calStfLis.rename_column("t", "t0")

    calStfLis = calStfLis["name", "kp", "xarc", "yarc", 
                      "xerr", "yerr", "vx", 
                      "vy", "vxerr", "vyerr",
                     "t0", "use", "r2d"]

    def replaceStarWithS00(lst):
        lnew = []
        for s in lst:
            if s.find("star") >= 0:
                s = s.replace("star_", "")
                s_new = "S"
                for i in range(3 - len(s)):
                    s_new += "0"
                lnew.append(s_new + s)
            else:
                 lnew.append(s)
        return lnew

    calStfLis["name"] = replaceStarWithS00(calStfLis["name"])
    #calStfLis.pprint()

    label_file = "../../source_list/" + target + "_label.dat"
    _out = open(label_file, 'w+')
    _out.write('#')
    ascii.write(calStfLis, _out, format='fixed_width', delimiter=' ', bookend=False)
    _out.close()


# Run align RMS, generates phot_error, pos_error, etc files.
def align():
    analysisObject = analysis.mb19284(epoch, filt)
    analysisObject.alignCombo()

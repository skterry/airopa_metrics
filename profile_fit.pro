;+
; :Description:
;     Perform photometry and astrometry on an image using AIROPA in "legacy",
;     "single PSF" and "variable PSF" modes.
;
; :Params:
;     img : in, required, type=STRING
;       Input image (without extension)
;     dir_test : in, required, type=STRING
;       Directory of the test
;     year : in, required, type=INT
;       Year of the phase map
;     filter_name : in, required, type=STRING
;       Name of the filter to be used
;     atm_inst : in, required, type=STRING
;       Atmospheric or instrument-only OTF
;     parang : in, required, type=FLOAT
;       Parallactic angle of the instrument (deg)
;     rotposn : in, required, type=FLOAT
;       Rotator position (deg)
;     pa : in, required, type=INT
;       Position angle
;     ref_star : in, required, type=STRING
;       Name of the PSF reference star
;
; :Keywords:
;     psf_const : in, required, type=INT, default=0
;       Use a constant PSF (0 or 1)
;     psf_load: in, required, type=INT, default=0
;       Load the PSF from file (0 or 1)
;     corr : in, optional, type=FLOAT, default=0.8
;       Correlation value
;     deblend : in, optional, type=INT, default=0
;       Deblending stars (0 or 1)
;     trimfake : in, optional, type=INT, default=1
;       Remove faint stars too close to others (0 or 1)
;     fixpsf : in, optional, type=INT, default=0
;       Clean the PSF (0 or 1)
;     back_clip : in, optional, type=INT, default=1
;       noise-dominated areas of the PSF are clipped to zero
;     weighted : in, optional, type=INT, default=1
;       Weight PSF stars (0 or 1)
;     backboxFWHM : in, optional, type=FLOAT, default=25.
;       Box size for the background
;     n_grid_over : in, optional, type=INT, default=3
;       Grid oversampling
;     part_size : in, optional, type=INT, default=102
;       Grid pitch (px)
;     flat : in, optional, type=INT, default=1
;       Flat background (0 or 1)
;     subtract : in, optional, type=INT, default=1
;       Subtract threshold (0 or 1)
;     dimm_file : in, required, type=STRING
;       DIMM file name
;     mass_file : in, required, type=STRING
;       MASS file name
;     atm_corr_file : in, required, type=STRING
;       File with the correction for resampled atmospheric OTF
;     on_sky : in, optional, type=INT, default=0
;        On-sky observation, meaning that the atmospheric OTF has not been
;        saved (0 or 1)
;     force : in, optional, type=INT, default=0
;        Force initial catalog (0 or 1)
;     x_force : in, optional, type=FLOAT Array[1, n]
;        Initial forced catalog x positions
;     y_force : in, optional, type=FLOAT Array[1, n]
;        Initial forced catalog y positions
;     f_force : in, optional, type=FLOAT Array[1, n]
;        Initial forced catalog fluxes
;     conf_file: in, required, type=STRING
;       AIROPA config file name
;
; :Uses:
;     bool_to_str
;     find_stf
;-
PRO profile_fit, img, dir_test, year, filter_name, atm_inst, parang, rotposn, $
  pa, ref_star, psf_const=psf_const, psf_load=psf_load, corr=corr, $
  deblend=deblend, trimfake=trimfake, fixpsf=fixpsf, back_clip=back_clip, $
  weighted=weighted, backboxFWHM=backboxFWHM, n_grid_over=n_grid_over, $
  part_size=part_size, flat=flat, subtract=subtract, dimm_file=dimm_file, $
  mass_file=mass_file, atm_corr_file=atm_corr_file, conf_file=conf_file, $
  on_sky=on_sky, force=force, x_force=x_force, y_force=y_force, f_force=f_force, $
  fix_Psf_Cos = fix_Psf_Cos, fix_Psf_Smooth = fix_Psf_Smooth, fix_Psf_HaloClip = fix_Psf_HaloClip, $
  fix_Psf_Trim = fix_Psf_Trim, fix_Psf_maskrad = fix_Psf_maskrad, fix_Psf_nSigmaStart = fix_Psf_nSigmaStart
  
  ; Set the keywords' default values
  IF N_ELEMENTS(psf_const) EQ 0 THEN (psf_const = 0)
  IF N_ELEMENTS(psf_load) EQ 0 THEN (psf_load = 0)
  IF NOT KEYWORD_SET(corr) THEN (corr = 0.8) ; changed from 0.6
  IF N_ELEMENTS(deblend) EQ 0 THEN (deblend = 0)
  IF N_ELEMENTS(trimfake) EQ 0 THEN (trimfake = 0)
  IF N_ELEMENTS(fixpsf) EQ 0 THEN (fixpsf = 1) ; changed from 1
  IF N_ELEMENTS(back_clip) EQ 0 THEN (back_clip = 0)
  IF N_ELEMENTS(weighted) EQ 0 THEN (weighted = 0)
  IF NOT KEYWORD_SET(backboxFWHM) THEN (backboxFWHM = 25)
  IF NOT KEYWORD_SET(n_grid_over) THEN (n_grid_over = 3)
  IF NOT KEYWORD_SET(part_size) THEN (part_size = 102)
  IF NOT KEYWORD_SET(clean_rad) THEN (clean_rad = 25)
  IF N_ELEMENTS(flat) EQ 0 THEN (flat = 1) ; changed from 1
  IF N_ELEMENTS(subtract) EQ 0 THEN (subtract = 1)
  IF N_ELEMENTS(on_sky) EQ 0 THEN (on_sky = 1) ; changed from 0
  IF N_ELEMENTS(force) EQ 0 THEN (force = 0)
  IF NOT KEYWORD_SET(atm_corr_file) THEN (atm_corr_file = $
    'ratio_atm_corr_k_101_to_301.fits')
  IF NOT KEYWORD_SET(conf_file) THEN (conf_file = 'airopa.config')
  IF NOT KEYWORD_SET(fix_Psf_Cos) THEN (fix_Psf_Cos = 0)
  IF NOT KEYWORD_SET(fix_Psf_Smooth) THEN (fix_Psf_Smooth = 0)
  IF NOT KEYWORD_SET(fix_Psf_HaloClip) THEN (fix_Psf_HaloClip = 0)
  IF NOT KEYWORD_SET(fix_Psf_Trim) THEN (fix_Psf_Trim = 0)
  IF NOT KEYWORD_SET(fix_Psf_maskrad) THEN (fix_Psf_maskrad = 0)
  IF NOT KEYWORD_SET(fix_Psf_nSigmaStart) THEN (fix_Psf_nSigmaStart = 0)
  
  ; Generate a header file
  CALDAT, SYSTIME(/JULIAN), month_cal, day_cal, year_cal
  IF STRCMP(atm_inst, 'inst') THEN phase_type = 'instrument-only'
  IF STRCMP(atm_inst, 'atm') THEN phase_type = 'atmospheric'
  head_txt = $
    'Image:                    ' + (img + '.fits') + STRING(10B) + $
    'Folder:                   ' + dir_test + STRING(10B) + $
    'Date analyzed:            ' + STRTRIM(STRING(year_cal, FORMAT='(I)'), 1) $
    + '-' + STRTRIM(STRING(month_cal, FORMAT='(I)'), 1) + '-' + $
    STRTRIM(STRING(day_cal, FORMAT='(I)'), 1) + STRING(10B) + STRING(10B) + $
    '--- PARAMETERS ---' + STRING(10B) + $
    'Phase map year:           ' + STRTRIM(STRING(year), 1) + STRING(10B) + $
    'Phase map type:           ' + phase_type + STRING(10B) + $
    'Filter:                   ' + filter_name + STRING(10B) + $
    'Parallactic angle (deg):  ' + STRTRIM(STRING(parang, FORMAT='(F5.1)'), 1) $
    + STRING(10B) + $
    'Rotator position (deg):   ' + $
    STRTRIM(STRING(rotposn, FORMAT='(F5.1)'), 1) + STRING(10B) + $
    'Position angle (deg):     ' + STRTRIM(STRING(pa), 1) + STRING(10B) + $
    'Reference star:           ' + ref_star + STRING(10B) + $
    'Constant PSF:             ' + bool_to_str(psf_const) + STRING(10B) + $
    'Load PSF:                 ' + bool_to_str(psf_load) + STRING(10B) + $
    'Correlation:              ' + STRTRIM(STRING(corr, FORMAT='(F0.2)'), 1) + $
    STRING(10B) + $
    'Deblend:                  ' + bool_to_str(deblend) + STRING(10B) + $
    'Trim fake:                ' + bool_to_str(trimfake) + STRING(10B) + $
    'Fix PSF:                  ' + bool_to_str(fixpsf) + STRING(10B) + $
    'Background clip:          ' + bool_to_str(back_clip) + STRING(10B) + $
    'Weighted:                 ' + bool_to_str(weighted) + STRING(10B) + $
    'Backbox FWHM (px):        ' + STRTRIM(STRING(backboxFWHM, $
    FORMAT='(F0.1)'), 1) + STRING(10B) + $
    'Clean radius (px):        ' + STRTRIM(STRING(clean_rad, FORMAT='(F0.1)'), $
    1) + STRING(10B) + $
    'Grid oversampling:        ' + STRTRIM(STRING(n_grid_over, FORMAT='(I)'), $
    1) + STRING(10B) + $
    'Grid pitch (px)  :        ' + STRTRIM(STRING(part_size, FORMAT='(I)'), 1) $
    + STRING(10B) + $
    'Flat background:          ' + bool_to_str(flat) + STRING(10B) + $
    'Subtract threshold:       ' + bool_to_str(subtract) + STRING(10B) + $
    'Force initial catalog:    ' + bool_to_str(force) + STRING(10B) + $
    'Atm OTF correction file:  ' + atm_corr_file + STRING(10B) + $
    'Configuration file:       ' + conf_file + STRING(10B) + $
    'DIMM file:                ' + dimm_file + STRING(10B) + $
    'MASS file:                ' + mass_file + STRING(10B)
  OPENW, 2, FILEPATH(('params_fit_' + img + '.txt'), ROOT_DIR=dir_test, $
    SUBDIR='fit')
  PRINTF, 2, head_txt
  CLOSE, 2
  
  ; Use different AIROPA modes
  psf_const_str = ''
  IF psf_const THEN psf_const_str = '_flat'
  mode_name = ['legacy', 'single', 'variable']
  mode_opt = [[1, 0], [0, 0], [0, 1]] ; Options [legacy, aoopt] for "legacy",
    ; "single PSF" and "variable PSF"
  IF psf_load THEN BEGIN
    makePsf = 0
    modes = [1, 1]
  ENDIF ELSE BEGIN
    makePsf = 1
    modes = [0, 2]
  ENDELSE
  IF STRCMP(atm_inst, 'inst') THEN BEGIN
    makeGrid_flag = 0
  ENDIF ELSE BEGIN
    makeGrid_flag = 1
  ENDELSE
  IF on_sky THEN BEGIN
    inst_grid_name_real_val = 0
    inst_grid_name_imag_val = 0
  ENDIF ELSE BEGIN
    inst_otf_grid_real = 'otf_ratio_grid_real' + psf_const_str + '.fits'
    inst_otf_grid_imag = 'otf_ratio_grid_imag' + psf_const_str + '.fits'
    inst_grid_name_real_val = FILEPATH(inst_otf_grid_real, $
      ROOT_DIR=GETENV('AIROPA_DATA_PATH'), SUBDIR=['otf_grids', $
      STRTRIM(STRING(year), 1), filter_name, 'inst', STRTRIM(STRING(pa), 1)])
    inst_grid_name_imag_val = FILEPATH(inst_otf_grid_imag, $
      ROOT_DIR=GETENV('AIROPA_DATA_PATH'), SUBDIR=['otf_grids', $
      STRTRIM(STRING(year), 1), filter_name, 'inst', STRTRIM(STRING(pa), 1)])
    otf_grid_name_real_val = FILEPATH(inst_otf_grid_real, $
      ROOT_DIR=GETENV('AIROPA_DATA_PATH'), SUBDIR=['otf_grids', $
      STRTRIM(STRING(year), 1), filter_name, 'inst', STRTRIM(STRING(pa), 1)])
    otf_grid_name_imag_val = FILEPATH(inst_otf_grid_imag, $
      ROOT_DIR=GETENV('AIROPA_DATA_PATH'), SUBDIR=['otf_grids', $
      STRTRIM(STRING(year), 1), filter_name, 'inst', STRTRIM(STRING(pa), 1)])
  ENDELSE
  gridposFile_value = FILEPATH('grid_pos.fits', $
    ROOT_DIR=GETENV('AIROPA_TEST_DATA'), SUBDIR=['psf_grids', $
    STRTRIM(STRING(year), 1), filter_name, 'inst', STRTRIM(STRING(pa), 1)])
  JOURNAL, FILEPATH(('journal_fit_' + img + '.txt'), ROOT_DIR=dir_test, $
    SUBDIR='fit')
  FOR i_mode = modes[0], modes[1] DO BEGIN;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;modes[0]
    ; Perform the profile fitting
    FILE_MKDIR, FILEPATH(PATH_SEP(), ROOT_DIR=dir_test, $
      SUBDIR=['fit', mode_name[i_mode]])
    CD, FILEPATH(PATH_SEP(), ROOT_DIR=dir_test, SUBDIR='image'), CURRENT=old_dir
    find_stf, $
      FILEPATH((img + '.fits'), ROOT_DIR=dir_test, SUBDIR='image'), $
      corr, $
      cooStar=ref_star, $
      gsStar='', $
      ttStar='TTstar', $
      starlist=FILEPATH('psf_stars.dat', ROOT_DIR=dir_test, SUBDIR='image'), $
      flat=flat, $
      subtract=subtract, $
      fixPsf=fixPsf, $
      trimfake=trimfake, $
      deblend=deblend, $
      back_clip=back_clip, $
      weighted=weighted, $
      backboxFWHM=backboxFWHM, $
      clean_rad=clean_rad, $
      n_grid_over=n_grid_over, $
      part_size=part_size, $
      makePsf=makePsf, $
      makeGrid=makeGrid_flag, $
      legacy=mode_opt[0, i_mode], $
      aoopt=mode_opt[1, i_mode], $
      inst_grid_name_real=inst_grid_name_real_val, $
      inst_grid_name_imag=inst_grid_name_imag_val, $
      otfgridFile_real=otf_grid_name_real_val, $
      otfgridFile_imag=otf_grid_name_imag_val, $
      posFile=gridposFile_value, $
      phase_map_folder=FILEPATH('', ROOT_DIR=GETENV('AIROPA_DATA'), $
      SUBDIR=['phase_maps', STRTRIM(STRING(year), 1)]), $
      dimm_file=FILEPATH(dimm_file, ROOT_DIR=dir_test, SUBDIR='image'), $
      mass_file=FILEPATH(mass_file, ROOT_DIR=dir_test, SUBDIR='image'), $
      corr_file=FILEPATH(atm_corr_file, ROOT_DIR=GETENV('AIROPA_DATA'), $
      SUBDIR='ref_files'), $
      config_file=FILEPATH(conf_file, ROOT_DIR=GETENV('AIROPA_DATA'), $
      SUBDIR='ref_files'), $
      force=force, $
      x_force=x_force, $
      y_force=y_force, $
      f_force=f_force, $
      makeRes=1, $
      makeStars=1, $
      debug=0, $
      save_otf=0, $
      fix_Psf_Cos = fix_Psf_Cos, $
      fix_Psf_Smooth = fix_Psf_Smooth, $
      fix_Psf_HaloClip = fix_Psf_HaloClip, $
      fix_Psf_Trim = fix_Psf_Trim, $
      fix_Psf_maskrad = fix_Psf_maskrad, $
      fix_Psf_nSigmaStart = fix_Psf_nSigmaStart
    CD, old_dir
    
    ; Move the output files to their relative output folder
    IF (i_mode EQ 2) THEN BEGIN
      psf_str = '_on_axis'
    ENDIF ELSE BEGIN
      psf_str = ''
    ENDELSE
    FILE_MOVE, FILEPATH([(img + '_res.fits'), (img + '_stars.fits')], $
      ROOT_DIR=dir_test, SUBDIR='image'), FILEPATH('', ROOT_DIR=dir_test, $
      SUBDIR=['fit', mode_name[i_mode]]), /OVERWRITE
    FILE_MOVE, FILEPATH([(img + '_' + STRMID(STRTRIM(corr, 1), 0, 3) + $
      '_stf.txt'), (img + '_' + STRMID(STRTRIM(corr, 1), 0, 3) + $
      '_metrics.txt'), (img + '_' + STRMID(STRTRIM(corr, 1), 0, 3) + $
      '_stf.lis')], ROOT_DIR=dir_test, SUBDIR='image'), FILEPATH('', $
      ROOT_DIR=dir_test, SUBDIR=['fit', mode_name[i_mode]]), /OVERWRITE
    IF NOT psf_load THEN BEGIN
      FILE_MOVE, FILEPATH([(img + '_back.fits'), (img + psf_str + $
        '_psf.fits')], ROOT_DIR=dir_test, SUBDIR='image'), FILEPATH('', $
        ROOT_DIR=dir_test, SUBDIR=['fit', mode_name[i_mode]]), /OVERWRITE
    ENDIF
    IF (i_mode EQ 2) THEN BEGIN
      FILE_MOVE, FILEPATH((img + '_psf_grid.fits'), ROOT_DIR=dir_test, $
        SUBDIR='image'), FILEPATH('', ROOT_DIR=dir_test, SUBDIR=['fit', $
        mode_name[i_mode]]), /OVERWRITE
    ENDIF
    IF (i_mode EQ 2) and (makeGrid_flag EQ 1) THEN BEGIN
      FILE_MOVE, FILEPATH((img + '_grid_pos.fits'), ROOT_DIR=dir_test, $
        SUBDIR='image'), FILEPATH('', ROOT_DIR=dir_test, SUBDIR=['fit', $
        mode_name[i_mode]]), /OVERWRITE
    ENDIF
  ENDFOR
  JOURNAL
END

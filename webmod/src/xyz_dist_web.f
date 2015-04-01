***********************************************************************
c     xyz_dist_web.f: Distribute temperatures and precip to MRUs using
c                      Lauren Hay's XYZ methodology.
c
c                      Converted by Steve Markstrom
c                      Wed Feb 10 15:16:04 MST 1999
c              revised Wed Mar 17 15:48:52 MST 1999
c              revised Mon Aug 30 16:47:07 MDT 1999
c              revised Mon Aug 30 16:47:07 MDT 1999
c              revised Wed Mar  8 09:06:18 MST 2000
c              revised Thu Feb  3 10:00:00 MST 2005
c                 (version 1.3 check by markstro -- no declpri needed)
c
c              revised and converted from xyz_dist_soltsta.f
c                 for use in webmod by Rick Webb Tues Mar 15 2005
c              Added version control
c
c***********************************************************************

C
C nMRU - number of mrus     (RW: webmod uses MRUs instead of HRUs)
C
C temp_nsta - number of temperature stations used
C temp_nuse (temp_nsta) - indicies of temp stations used
C ntemp - total number of temp stations
C
C rain_nsta - number of precip stations used
C rain_nuse (rain_nsta) - indicies of rain stations used
C nrain - total number of precip stations

c***********************************************************************
c 
c     xyzdecl - set up parameters for temperature computations
c

      integer function xyzdecl (tmax_f, tmin_f, temp_f,
     +                           tmax_c, tmin_c, temp_c,
     +                           basin_tmax_c, basin_tmin_c,
     +                           basin_tmax_f, basin_tmin_f,
     $                           basin_temp_c, basin_temp_f,
     +                           mru_ppt, mru_dep,
     $                           basin_ppt, basin_dep,
     $                           newsnow,
     +                           mru_rain, mru_snow, pptmix,
     +                           prmx, tmax_rain_sta, tmin_rain_sta,
     +                           is_rain_day,solrad_tmax)

      include 'fmodules.inc'

      integer newsnow(MAXMRU), pptmix(MAXMRU)
      integer is_rain_day
      real tmax_f(MAXMRU), tmin_f(MAXMRU), temp_f(MAXMRU)
      real tmax_c(MAXMRU), tmin_c(MAXMRU), temp_c(MAXMRU)
      real mru_ppt(MAXMRU),mru_dep(MAXMRU)
      real basin_tmax_c, basin_tmin_c, basin_temp_c
      real basin_tmax_f, basin_tmin_f, basin_temp_f
      real basin_ppt,basin_dep
      real mru_snow(MAXMRU), mru_rain(MAXMRU)
      real prmx(MAXMRU)
      real tmax_rain_sta(MAXRAIN), tmin_rain_sta(MAXRAIN)
      real solrad_tmax

      xyzdecl = 1


C23456789112345678921234567893123456789412345678951234567896123456789712
      if(declvar('xyz_dist', 'solrad_tmax', 'one', 1, 'real',
     +     'Basin max temp adjusted to elevation of solrad station',
     +     'degrees F',
     +   solrad_tmax).ne.0) return

      if(declvar('xyz_dist', 'tmax_f', 'nmru', MAXMRU, 'real',
     +     'MRU adjusted daily maximum temperature',
     +     'degrees F',
     +   tmax_f).ne.0) return

      if(declvar('xyz_dist', 'is_rain_day', 'one', 1, 'integer',
     +     'Is it raining in the basin',
     +     'none',
     +   is_rain_day).ne.0) return

      if(declvar('xyz_dist', 'tmin_f', 'nmru', MAXMRU, 'real',
     +     'MRU adjusted daily minimum temperature',
     +     'degrees F',
     +   tmin_f).ne.0) return

      if(declvar('xyz_dist', 'temp_f', 'nmru', MAXMRU, 'real',
     +     'MRU adjusted daily average temperature',
     +     'degrees F',
     +   temp_f).ne.0) return


C      if(decl var('temp', 'tempf', 'nmru', MAXMRU, 'real',
C     +     'MRU adjusted temperature for timestep',
C     +     'degrees F',
C     +   tempf).ne.0) return

      if(declvar('xyz_dist', 'tmax_c', 'nmru', MAXMRU, 'real',
     +     'MRU adjusted daily maximum temperature',
     +     'degrees C',
     +   tmax_c).ne.0) return

      if(declvar('xyz_dist', 'tmin_c', 'nmru', MAXMRU, 'real',
     +     'MRU adjusted daily minimum temperature',
     +     'degrees C',
     +   tmin_c).ne.0) return

      if(declvar('xyz_dist', 'temp_c', 'nmru', MAXMRU, 'real',
     +     'MRU adjusted daily average temperature',
     +     'degrees F',
     +   temp_c).ne.0) return

C      if(decl var('temp', 'tempc', 'nmru', MAXMRU, 'real',
C     +     'MRU adjusted temperature for timestep',
C     +     'degrees C',
C     +   tempc).ne.0) return

      if(declvar('xyz_dist', 'basin_tmax_c', 'one', 1, 'real',
     +     'Basin area-weighted daily maximum temperature',
     +     'degrees C',
     +   basin_tmax_c).ne.0) return

      if(declvar('xyz_dist', 'basin_tmin_c', 'one', 1, 'real',
     +     'Basin area-weighted daily minimum temperature',
     +     'degrees C',
     +   basin_tmin_c).ne.0) return

      if(declvar('xyz_dist', 'basin_temp_c', 'one', 1, 'real',
     +     'Basin area-weighted temperature for timestep',
     +     'degrees C',
     +   basin_temp_c).ne.0) return

      if(declvar('xyz_dist', 'basin_tmax_f', 'one', 1, 'real',
     +     'Basin area-weighted daily maximum temperature',
     +     'degrees F',
     +   basin_tmax_f).ne.0) return

      if(declvar('xyz_dist', 'basin_tmin_f', 'one', 1, 'real',
     +     'Basin area-weighted daily minimum temperature',
     +     'degrees F',
     +   basin_tmin_f).ne.0) return

      if(declvar('temp', 'basin_temp_f', 'one', 1, 'real',
     +     'Basin area-weighted temperature for timestep',
     +     'degrees F',
     +   basin_temp_f).ne.0) return

      if(declvar('xyz_dist', 'mru_ppt', 'nmru', MAXMRU, 'real',
     +     'Adjusted precip on each MRU',
     +     'inches',
     +   mru_ppt).ne.0) return

      if(declvar('xyz_dist', 'mru_dep', 'nmru', MAXMRU, 'real',
     +     'Adjusted precip+irrig on each MRU',
     +     'inches',
     +   mru_dep).ne.0) return

C      if(decl var('temp', 'basin_temp', 'one', 1, 'real',
C     +     'Basin area-weighted temperature for timestep',
C     +     'degrees',
C     +   basin_temp).ne.0) return

      if(declvar('xyz_dist', 'basin_ppt', 'one', 1, 'real',
     +     'Area weighted adjusted average precip for basin',
     +     'inches',
     +   basin_ppt).ne.0) return

      if(declvar('xyz_dist', 'basin_dep', 'one', 1, 'real',
     +     'Area-weighted adjusted average precip+irrig for basin',
     +     'inches',
     +   basin_dep).ne.0) return

C DANGER - Not sure what to do about this one.  For right now
C          I'm setting basin_ppt and basin_obs_ppt to the same
C          variable.  In the precip_prms module, basin_obs_ppt
C          seems to be the area weighted precip average before
C          the correction factor is applied.  In this module,
C          the correction "error" is applied to the station
C          precip rather than the mru precip.
c Eliminated, This is basin_ppt or basin_dep - RW
c
c$$$      if(decl var('xyz_dist', 'basin_obs_ppt', 'one', 1, 'real',
c$$$     +     'Area weighted measured average precip for basin',
c$$$     +     'inches',
c$$$     +   basin_ppt).ne.0) return
c$$$
      if(declvar('xyz_dist', 'newsnow', 'nmru', MAXMRU, 'integer',
     +     'New snow on MRU, 0=no, 1=yes',
     +     'none',
     +   newsnow).ne.0) return

      if(declvar('precip', 'mru_rain', 'nmru', MAXMRU, 'real',
     +     'Computed rain on each MRU',
     +     'inches',
     +   mru_rain).ne.0) return

      if(declvar('precip', 'mru_snow', 'nmru', MAXMRU, 'real',
     +     'Computed snow on each MRU',
     +     'inches',
     +   mru_snow).ne.0) return

      if(declvar('precip', 'pptmix', 'nmru', MAXMRU, 'integer',
     +     'Precip mixture - 0=no, 1=yes',
     +     'none',
     +   pptmix).ne.0) return

      if(declvar('precip', 'prmx', 'nmru', MAXMRU, 'real',
     +     'Proportion of rain in a mixed event',
     +     'decimal percent; equal to -0.1 on days with no rain',
     +   prmx).ne.0) return

C23456789112345678921234567893123456789412345678951234567896123456789712
      if(declvar('xyz_dist', 'tmax_rain_sta', 'nrain', MAXRAIN, 'real',
     +     'Maximum temperature distributed to the precip stations',
     +     'degrees F',
     +   tmax_rain_sta).ne.0) return

      if(declvar('xyz_dist', 'tmin_rain_sta', 'nrain', MAXRAIN, 'real',
     +     'Minimum temperature distributed to the precip stations',
     +     'degrees F',
     +   tmin_rain_sta).ne.0) return

      if (declparam ('xyz_dist', 'mru_elev', 'nmru', 'real',
     +   '0.', '-300.', '10000',
     +   'Mean elevation for each MRU',
     +   'Mean elevation for each MRU',
     +   'meters').ne.0) return

      if (declparam ('xyz_dist', 'mru_x', 'nmru', 'real',
     +   '0.', '-10000000.', '10000000.',
     +   'X for each MRU (albers)',
     +   'X for each MRU (albers)',
     +   'meters').ne.0) return

      if (declparam ('xyz_dist', 'mru_y', 'nmru', 'real',
     +   '0.', '-10000000.', '10000000.',
     +   'Y for each MRU (albers)',
     +   'Y for each MRU (albers)',
     +   'meters').ne.0) return

      if (declparam ('xyz_dist', 'max_lapse', 'nlapse,nmonths',
     +   'real',
     +   '0.', '-10.', '10',
     +   'Maximum temperature lapse rate',
     +   'Maximum temperature lapse rate',
     +   'degree celsius per meter').ne.0) return

      if (declparam ('xyz_dist', 'min_lapse', 'nlapse,nmonths',
     +   'real',
     +   '0.', '-10.', '10',
     +   'Minimum temperature lapse rate',
     +   'Minimum temperature lapse rate',
     +   'degree celsius per meter').ne.0) return

      if (declparam ('xyz_dist', 'ppt_lapse', 'nlapse,nmonths',
     +   'real',
     +   '0.', '-10.', '10',
     +   'Precipitation lapse rate',
     +   'Precipitation lapse rate',
     +   'inch per meter').ne.0) return

      if (declparam ('xyz_dist', 'tsta_elev', 'ntemp', 'real',
     +   '0', '-300.', '10000.',
     +   'Temperature station elevation',
     +   'Elevation of each temperature '//
     +   'measurement station',
     +   'meters').ne.0) return

      if (declparam ('xyz_dist', 'tsta_x', 'ntemp', 'real',
     +   '0.', '-10000000.', '10000000.',
     +   'X for each temperature station (albers)',
     +   'X for each temperature station (albers)',
     +   'meters').ne.0) return

      if (declparam ('xyz_dist', 'tsta_y', 'ntemp', 'real',
     +   '0.', '-10000000.', '10000000.',
     +   'Y for each temperature station (albers)',
     +   'Y for each temperature station (albers)',
     +   'meters').ne.0) return

      if (declparam ('xyz_dist', 'psta_elev', 'nrain', 'real',
     +   '0', '-300.', '30000.',
     +   'Precip station elevation',
     +   'Elevation of each precip '//
     +   'measurement station',
     +   'meters').ne.0) return

      if (declparam ('xyz_dist', 'psta_x', 'nrain', 'real',
     +   '0.', '-10000000.', '10000000.',
     +   'X for each precip station (albers)',
     +   'X for each precip station (albers)',
     +   'meters').ne.0) return

      if (declparam ('xyz_dist', 'psta_y', 'nrain', 'real',
     +   '0.', '-10000000.', '10000000.',
     +   'Y for each precip station (albers)',
     +   'Y for each precip station (albers)',
     +   'meters').ne.0) return

      if (declparam ('xyz_dist', 'tsta_nuse', 'ntemp', 'integer',
     +   '1', '0', '1',
     +   '0 = station not used  1 = station used',
     +   '0 = station not used  1 = station used',
     +   'none').ne.0) return

C23456789112345678921234567893123456789412345678951234567896123456789712
      if (declparam ('xyz_dist', 'solrad_elev', 'one', 'real',
     +   '1000.0', '0.0', '10000.0',
     +   'Elevation of the solrad station used for the DD curves.',
     +   'Elevation of the solrad station used for the DD curves.',
     +   'meters').ne.0) return

C23456789112345678921234567893123456789412345678951234567896123456789712
      if (declparam ('xyz_dist', 'psta_nuse', 'nrain', 'integer',
     +   '1', '0', '1',
     +   'The subset of precip stations used in the distribution'//
     +   ' regression.  0 = station not used  1 = station used',
     +   'The subset of precip stations used in the distribution'//
     +   ' regression.  0 = station not used  1 = station used',
     +   'none').ne.0) return

C23456789112345678921234567893123456789412345678951234567896123456789712
      if (declparam ('xyz_dist', 'psta_freq_nuse', 'nrain', 'integer',
     +   '1', '0', '1',
     +   'The subset of precip stations used to determine if '//
     +   ' there is precip in the basin.  0 = station not used '//
     +   ' 1 = station used',
     +   'The subset of precip stations used to determine if '//
     +   ' there is precip in the basin.  0 = station not used '//
     +   ' 1 = station used',
     +   'none').ne.0) return

      if(declparam('xyz_dist', 'mru_area', 'nmru', 'real',
     +   '1.0', '0.01', '1e+09',
     +   'MRU area',
     +   'MRU area',
     +   'km2').ne.0) return

      if(declparam('xyz_dist', 'basin_area', 'one', 'real',
     +   '1.0', '0.01', '1e+09',
     +   'Total basin area',
     +   'Total basin area',
     +   'km2').ne.0) return

      if (declparam ('xyz_dist', 'tsta_month_max',
     +   'ntemp,nmonths', 'real',
     +   '0.', '-100.', '200.',
     +   'Average monthly maximum temp at each station',
     +   'Average monthly maximum temp at each station',
     +   'degrees').ne.0) return

      if (declparam ('xyz_dist', 'tsta_month_min',
     +   'ntemp,nmonths', 'real',
     +   '0.', '-100.', '200.',
     +   'Average monthly minimum temp at each station',
     +   'Average monthly minimum temp at each station',
     +   'degrees').ne.0) return

      if (declparam ('xyz_dist', 'psta_month_ppt',
     +   'nrain,nmonths', 'real',
     +   '0.', '0.', '200.',
     +   'Average monthly precip at each station',
     +   'Average monthly precip at each station',
     +   'inches').ne.0) return

      if (declparam ('xyz_dist', 'adjust_snow',
     +   'nmonths', 'real',
     +   '0.01', '0.0', '1.0',
     +   'Downscaling % adjustment for snow',
     +   'Downscaling % adjustment for snow',
     +   'decimal percent').ne.0) return

      if (declparam ('xyz_dist', 'adjust_rain',
     +   'nmonths', 'real',
     +   '0.01', '0.0', '1.0',
     +   'Downscaling % adjustment for rain',
     +   'Downscaling % adjustment for rain',
     +   'decimal percent').ne.0) return

      if(declparam('xyz_dist', 'tmax_allrain_c', 'nmonths', 'real',
     +   '5.', '0.', '15.',
     +   'Precip all rain if tmax_c above this value',
     +   'If MRU maximum temperature exceeds this value, '//
     +   'precipitation is assumed to be rain.','degrees celsius')
     +   .ne.0) return

      if(declparam('xyz_dist', 'rain_code', 'nmonths', 'integer',
     +   '2', '1', '5',
     +   'Code indicating rule for precip station use',
     +   'Code indicating rule for precip station use: '//
     +   '1 = only precip if the regression stations have precip,'//
     +   '2 = only precip if any station in the basin has precip,'//
     +   '3 = precip if xyz says so,'//
     +   '4 = only precip if rain_day variable is set to 1,'//
     +   '5 = only precip if psta_freq_nuse stations see precip',
     +   'none')
     +   .ne.0) return

      if(declparam('xyz_dist', 'tmax_allsnow_c', 'one', 'real',
     +   '0', '-10.', '10.',
     +   'All snow if tmax_c< this value; all rain if tmin_c>this '//
     $   'value.',' If MRU maximum temperature is below this value, '//
     +   'precipitation is assumed to be snow; alternately, if MRU '//
     +   'minimum temperature is above this value, precipitation is '//
     +   'assumed to be all rain.','degrees celsius').ne.0) return

      if(declparam('xyz_dist', 'adjmix_rain', 'nmonths', 'real',
     +   '1.', '0.', '3.',
     +   'Adjustment factor for rain in a rain/snow mix',
     +   'Monthly factor to adjust rain proportion in a mixed '//
     +   'rain/snow event',
     +   'none').ne.0) return

      if(declparam('xyz_dist', 'x_add', 'one', 'real',
     +   '0.0', '-10000000.0', '10000000.0',
     +   'X additive term for climate station transform',
     +   'X additive term for climate station transform',
     +   'meters')
     +   .ne.0) return

      if(declparam('xyz_dist', 'x_div', 'one', 'real',
     +   '0.0', '-10000000.0', '10000000.0',
     +   'X divisor term for climate station transform',
     +   'X divisor term for climate station transform',
     +   'meters')
     +   .ne.0) return

      if(declparam('xyz_dist', 'y_add', 'one', 'real',
     +   '0.0', '-10000000.0', '10000000.0',
     +   'Y additive term for climate station transform',
     +   'Y additive term for climate station transform',
     +   'meters')
     +   .ne.0) return

      if(declparam('xyz_dist', 'y_div', 'one', 'real',
     +   '0.0', '-10000000.0', '10000000.0',
     +   'Y divisor term for climate station transform',
     +   'Y divisor term for climate station transform',
     +   'meters')
     +   .ne.0) return

      if(declparam('xyz_dist', 'z_add', 'one', 'real',
     +   '0.0', '-10000000.0', '10000000.0',
     +   'Z additive term for climate station transform',
     +   'Z additive term for climate station transform',
     +   'meters')
     +   .ne.0) return

      if(declparam('xyz_dist', 'z_div', 'one', 'real',
     +   '0.0', '-10000000.0', '10000000.0',
     +   'Z divisor term for climate station transform',
     +   'Z divisor term for climate station transform',
     +   'meters')
     +   .ne.0) return

      if(declparam('xyz_dist', 'tmax_add', 'one', 'real',
     +   '0.0', '-10000000.0', '10000000.0',
     +   'Max temp additive term for climate station transform',
     +   'Max temp additive term for climate station transform',
     +   'degrees C')
     +   .ne.0) return

      if(declparam('xyz_dist', 'tmax_div', 'one', 'real',
     +   '0.0', '-10000000.0', '10000000.0',
     +   'Max temp divisor term for climate station transform',
     +   'Max temp divisor term for climate station transform',
     +   'degrees C')
     +   .ne.0) return

      if(declparam('xyz_dist', 'tmin_add', 'one', 'real',
     +   '0.0', '-10000000.0', '10000000.0',
     +   'Min temp additive term for climate station transform',
     +   'Min temp additive term for climate station transform',
     +   'degrees C')
     +   .ne.0) return

      if(declparam('xyz_dist', 'tmin_div', 'one', 'real',
     +   '0.0', '-10000000.0', '10000000.0',
     +   'Min temp divisor term for climate station transform',
     +   'Min temp divisor term for climate station transform',
     +   'degrees C')
     +   .ne.0) return

      if(declparam('xyz_dist', 'ppt_add', 'one', 'real',
     +   '0.0', '-10000000.0', '10000000.0',
     +   'Precip additive term for climate station transform',
     +   'Precip additive term for climate station transform',
     +   'inches')
     +   .ne.0) return

      if(declparam('xyz_dist', 'ppt_div', 'one', 'real',
     +   '0.0', '-10000000.0', '10000000.0',
     +   'Precip divisor term for climate station transform',
     +   'Precip divisor term for climate station transform',
     +   'inches')
     +   .ne.0) return

      if(declparam('temp', 'tmax_adj', 'nmru', 'real',
     +   '0.0', '-10.', '10.0',
     +   'MRU maximum temperature adjustment',
     +   'Adjustment, in deg C, '//
     +   'to MRU maximum temperature '//
     +   'based on slope and aspect of MRU',
     +   'degrees C').ne.0) return

      if(declparam('temp', 'tmin_adj', 'nmru', 'real',
     +   '0.0', '-10.0', '10.0',
     +   'MRU minimum temperature adjustment',
     +   'Adjustment, in deg C '//
     +   'to MRU minimum temperature '//
     +   'based on slope and aspect of MRU',
     +   'degrees C').ne.0) return


      xyzdecl = 0

      return
      end

c***********************************************************************
c
c     xyzinit - Initialize xyz_dist module - get parameter values,
c

C23456789112345678921234567893123456789412345678951234567896123456789712
      integer function xyzinit (nMRU, MRUx, MRUy, MRUelev, max_lapse,
     +                           min_lapse, temp_nsta, rain_STAx,
     +                           rain_STAy, rain_STAelev,
     +                           temp_STAx,
     +                           temp_STAy, temp_STAelev,
     +                           mru_area, basin_area,
     +                           tmaxMTH, tminMTH,
     +                           adjust_snow, meantmax, meantmin,
     +                           rain_meanx, rain_meany, rain_meanz,
     +                           temp_meanx, temp_meany, temp_meanz,
     +                           ppt_lapse, meanppt,
     +                           tmax_allrain_c, tmax_allsnow_c,
     +                           pptMTH, adjmix_rain,
     +                           tmin_add, tmin_div, tmax_add,
     +                           tmax_div, ppt_add, ppt_div,
     +                           tmin_adj, tmax_adj,
     +                           rain_code, rain_nsta, ntemp, nrain,
     +                           temp_nuse, rain_nuse, psta_freq_nuse,
     + solrad_elev, basin_centroid_x, basin_centroid_y, adjust_rain,
     + nform)

      include 'fmodules.inc'

      integer nMRU, temp_nsta, rain_nsta
      integer ntemp, nrain, nform
      real adjust_snow(MAXMO), adjust_rain(MAXMO)
      integer tsta_nuse(MAXTEMP), psta_nuse(MAXRAIN)
      integer psta_freq_nuse(MAXRAIN)
      integer temp_nuse(MAXTEMP), rain_nuse(MAXRAIN)
      integer rain_code(MAXMO)
      real MRUx(MAXMRU), MRUy(MAXMRU), MRUelev(MAXMRU)
      real max_lapse(MAXLAPSE,MAXMO), min_lapse(MAXLAPSE,MAXMO)
      real ppt_lapse(MAXLAPSE,MAXMO)
      real temp_STAx(MAXTEMP),temp_STAy(MAXTEMP), temp_STAelev(MAXTEMP)
      real rain_STAx(MAXTEMP),rain_STAy(MAXTEMP), rain_STAelev(MAXTEMP)
      real mru_area(MAXMRU), basin_area
      real tmaxMTH(MAXTEMP,MAXMO), tminMTH(MAXTEMP,MAXMO)
      real pptMTH(MAXRAIN,MAXMO)
      real meantmax(MAXMO), meantmin(MAXMO)
      real rain_meany(MAXMO), rain_meanz(MAXMO), rain_meanx(MAXMO)
      real temp_meany(MAXMO), temp_meanz(MAXMO), temp_meanx(MAXMO)
      real meanppt(MAXMO)
      real tmax_allrain_c(MAXMO), tmax_allsnow_c
      real adjmix_rain(MAXMO)
      real x_add, x_div, y_add, y_div, z_add, z_div
      real tmin_add, tmin_div, tmax_add, tmax_div
      real ppt_add, ppt_div
      real tmin_adj(MAXMRU), tmax_adj(MAXMRU)
      real basin
      real solrad_elev, basin_centroid_x, basin_centroid_y

      integer i, j, m

      common /lauren/ x_add, x_div, y_add, y_div, z_add, z_div

      xyzinit = 1

      nMRU = getdim ('nmru')
      if(nMRU.eq.-1) return

      ntemp = getdim ('ntemp')
      nrain = getdim ('nrain')

      nform = getdim ('nform')
      if(nform.eq.-1) return

      if (getparam ('xyz_dist', 'solrad_elev', 1, 'real', solrad_elev)
     +   .ne.0) return

      if (getparam ('xyz_dist', 'mru_x', MAXMRU, 'real', MRUx)
     +   .ne.0) return

      if (getparam ('xyz_dist', 'mru_y', MAXMRU, 'real', MRUy)
     +   .ne.0) return

C23456789112345678921234567893123456789412345678951234567896123456789712

      if (getparam ('xyz_dist', 'mru_elev', MAXMRU, 'real',
     +    MRUelev) .ne.0) return

      if (getparam ('xyz_dist', 'max_lapse', MAXLAPSE*MAXMO,
     +    'real', max_lapse) .ne.0) return

      if (getparam ('xyz_dist', 'min_lapse', MAXLAPSE*MAXMO,
     +    'real', min_lapse) .ne.0) return

      if (getparam ('xyz_dist', 'ppt_lapse', MAXLAPSE*MAXMO,
     +    'real', ppt_lapse) .ne.0) return

      if (getparam ('xyz_dist', 'tsta_x', MAXTEMP,
     +    'real', temp_STAx) .ne.0) return

      if (getparam ('xyz_dist', 'tsta_y', MAXTEMP,
     +    'real', temp_STAy) .ne.0) return

      if (getparam ('xyz_dist', 'tsta_elev', MAXTEMP,
     +    'real', temp_STAelev) .ne.0) return

      if (getparam ('xyz_dist', 'psta_x', MAXRAIN,
     +    'real', rain_STAx) .ne.0) return

      if (getparam ('xyz_dist', 'psta_y', MAXRAIN,
     +    'real', rain_STAy) .ne.0) return

      if (getparam ('xyz_dist', 'psta_elev', MAXRAIN,
     +    'real', rain_STAelev) .ne.0) return

      if (getparam ('xyz_dist', 'tsta_nuse', MAXTEMP,
     +    'integer', tsta_nuse) .ne.0) return

      if (getparam ('xyz_dist', 'psta_nuse', MAXRAIN,
     +    'integer', psta_nuse) .ne.0) return

      if (getparam ('xyz_dist', 'psta_freq_nuse', MAXRAIN,
     +    'integer', psta_freq_nuse) .ne.0) return

      if(getparam('xyz_dist', 'mru_area', MAXMRU, 'real',
     +   mru_area).ne.0) return

      if(getparam('xyz_dist', 'basin_area', 1, 'real',
     +   basin_area).ne.0) return

      if (getparam ('xyz_dist', 'tsta_month_min',
     +    MAXTEMP*MAXMO, 'real', tminMTH) .ne.0) return

      if (getparam ('xyz_dist', 'tsta_month_max',
     +    MAXTEMP*MAXMO, 'real', tmaxMTH) .ne.0) return

      if (getparam ('xyz_dist', 'psta_month_ppt',
     +    MAXRAIN*MAXMO, 'real', pptMTH) .ne.0) return

      if (getparam ('xyz_dist', 'adjust_snow',
     +    MAXMO, 'real', adjust_snow) .ne.0) return

      if (getparam ('xyz_dist', 'adjust_rain',
     +    MAXMO, 'real', adjust_rain) .ne.0) return

      if(getparam('precip', 'tmax_allrain_c', MAXMO, 'real',
     +   tmax_allrain_c).ne.0) return

      if(getparam('precip', 'tmax_allsnow_c', 1, 'real',
     $     tmax_allsnow_c).ne.0) return

      if(getparam('precip', 'adjmix_rain', MAXMO, 'real', adjmix_rain)
     +   .ne.0) return

      if(getparam('xyz_dist', 'z_add', 1, 'real', z_add)
     +   .ne.0) return

      if(getparam('xyz_dist', 'z_div', 1, 'real', z_div)
     +   .ne.0) return

      if(getparam('xyz_dist', 'x_add', 1, 'real', x_add)
     +   .ne.0) return

      if(getparam('xyz_dist', 'x_div', 1, 'real', x_div)
     +   .ne.0) return

      if(getparam('xyz_dist', 'y_add', 1, 'real', y_add)
     +   .ne.0) return

      if(getparam('xyz_dist', 'y_div', 1, 'real', y_div)
     +   .ne.0) return

      if(getparam('xyz_dist', 'tmax_add', 1, 'real', tmax_add)
     +   .ne.0) return

      if(getparam('xyz_dist', 'tmax_div', 1, 'real', tmax_div)
     +   .ne.0) return

      if(getparam('xyz_dist', 'tmin_add', 1, 'real', tmin_add)
     +   .ne.0) return

      if(getparam('xyz_dist', 'tmin_div', 1, 'real', tmin_div)
     +   .ne.0) return

      if(getparam('xyz_dist', 'ppt_add', 1, 'real', ppt_add)
     +   .ne.0) return

      if(getparam('xyz_dist', 'ppt_div', 1, 'real', ppt_div)
     +   .ne.0) return

      if(getparam('xyz_dist', 'tmax_adj', MAXMRU, 'real',
     +    tmax_adj)
     +   .ne.0) return

      if(getparam('xyz_dist', 'tmin_adj', MAXMRU, 'real',
     +   tmin_adj)
     +   .ne.0) return

      if(getparam('xyz_dist', 'rain_code', MAXMO, 'integer',
     +   rain_code)
     +   .ne.0) return

c
c Compute basin centroid
c
C23456789112345678921234567893123456789412345678951234567896123456789712
      basin_centroid_x = 0.0
      basin_centroid_y = 0.0
      do i = 1, nMRU
         basin_centroid_x = basin_centroid_x + (mru_area(i) * MRUx(i))
         basin_centroid_y = basin_centroid_y + (mru_area(i) * MRUy(i))
      end do

      basin_centroid_x = basin_centroid_x / basin_area
      basin_centroid_y = basin_centroid_y / basin_area

c
c convert mru elevation from feet to meters and transform X and Y
c
      basin = 0.0
      do i = 1, nMRU
C         MRUelev(i) = MRUelev(i) * 0.3048000
         MRUelev(i) = (MRUelev(i) + z_add) / z_div
         MRUx(i) = (MRUx(i) + x_add) / x_div
         MRUy(i) = (MRUy(i) + y_add) / y_div
         basin = mru_area(i) + basin
      end do


c
c convert temperature station elevation from feet to meters and
c   transform X and Y
c
      temp_nsta = 0
      do i = 1, ntemp
C         STAelev(i) = STAelev(i) * 0.3048000
c        temp_STAelev(i) = (temp_STAelev(i) + z_add) / z_div
c        temp_STAx(i) = (temp_STAx(i) + x_add) / x_div
c        temp_STAy(i) = (temp_STAy(i) + y_add) / y_div
         if (tsta_nuse(i) .eq. 1) then
            temp_nsta = temp_nsta + 1
            temp_nuse(temp_nsta)=i
         end if 
      end do

      rain_nsta = 0
      do i = 1, nrain
C         STAelev(i) = STAelev(i) * 0.3048000
c        rain_STAelev(i) = (rain_STAelev(i) + z_add) / z_div
c        rain_STAx(i) = (rain_STAx(i) + x_add) / x_div
c        rain_STAy(i) = (rain_STAy(i) + y_add) / y_div
         if (psta_nuse(i) .eq. 1) then
            rain_nsta = rain_nsta + 1
            rain_nuse(rain_nsta)=i
         end if 
      end do

c
c calculate the station mean by month
c
C23456789112345678921234567893123456789412345678951234567896123456789712
      do m=1,12
         meanppt(m)=0.0
         meantmax(m)=0.0
         meantmin(m)=0.0
         rain_meanx(m)=0.0
         rain_meany(m)=0.0
         rain_meanz(m)=0.0
         temp_meanx(m)=0.0
         temp_meany(m)=0.0
         temp_meanz(m)=0.0

         do j=1,rain_nsta
            i = rain_nuse(j)

c           pptMTH(i,m) = (pptMTH(i,m) + ppt_add) / ppt_div

            meanppt(m)=meanppt(m)+pptMTH(i,m)
            rain_meanx(m)=rain_meanx(m)+rain_STAx(i)
            rain_meany(m)=rain_meany(m)+rain_STAy(i)
            rain_meanz(m)=rain_meanz(m)+rain_STAelev(i)
         end do

         do j=1,temp_nsta
            i = temp_nuse(j)

c           tminMTH(i,m) = (tminMTH(i,m) + tmin_add) / tmin_div
c           tmaxMTH(i,m) = (tmaxMTH(i,m) + tmax_add) / tmax_div

            meantmin(m)=meantmin(m)+tminMTH(i,m)
            meantmax(m)=meantmax(m)+tmaxMTH(i,m)
            temp_meanx(m)=temp_meanx(m)+temp_STAx(i)
            temp_meany(m)=temp_meany(m)+temp_STAy(i)
            temp_meanz(m)=temp_meanz(m)+temp_STAelev(i)
         end do

         meanppt(m)=meanppt(m)/float(rain_nsta)
         meantmin(m)=meantmin(m)/float(temp_nsta)
         meantmax(m)=meantmax(m)/float(temp_nsta)

         rain_meanx(m)=(rain_meanx(m)/float(rain_nsta))
         rain_meany(m)=(rain_meany(m)/float(rain_nsta))
         rain_meanz(m)=rain_meanz(m)/float(rain_nsta)

         temp_meanx(m)=(temp_meanx(m)/float(temp_nsta))
         temp_meany(m)=(temp_meany(m)/float(temp_nsta))
         temp_meanz(m)=temp_meanz(m)/float(temp_nsta)

      meanppt(m)=(meanppt(m) + ppt_add) / ppt_div
      meantmin(m)=(meantmin(m) + tmin_add) / tmin_div
      meantmax(m)=(meantmax(m)+ tmax_add) / tmax_div

      rain_meanx(m)=(rain_meanx(m) + x_add) / x_div
        rain_meany(m)=(rain_meany(m) + y_add) / y_div
      rain_meanz(m)=(rain_meanz(m) + z_add) / z_div

        temp_meanx(m)=(temp_meanx(m) + x_add) / x_div
        temp_meany(m)=(temp_meany(m) + y_add) / y_div
        temp_meanz(m)=(temp_meanz(m) + z_add) / z_div

C      write (*,*) m, temp_meanx(m), temp_meany(m), temp_meanz(m)
C      write (*,*) m, rain_meanx(m), rain_meany(m), rain_meanz(m)

C      write(*,*) max_lapse(1,m), max_lapse(2,m), max_lapse(3,m)
C      write(*,*) min_lapse(1,m), min_lapse(2,m), min_lapse(3,m)
C      write(*,*) ppt_lapse(1,m), ppt_lapse(2,m), ppt_lapse(3,m)

      end do

      xyzinit = 0

      return
      end

c***********************************************************************
c
c     xyzrun - Temperature calculation
c               calculates daily max and min temperature
c               using data from available stations
c               Outputs a daily max and min Temperature by mru elevation
c

C23456789112345678921234567893123456789412345678951234567896123456789712
      integer function xyzrun (nMRU, MRUx, MRUy, MRUelev, max_lapse,
     +                          min_lapse, temp_nsta,
     +                          rain_STAx, rain_STAy, rain_STAelev,
     +                          temp_STAx, temp_STAy, temp_STAelev,
     +                          tmax_f, tmin_f, temp_f,
     +                          tmax_c, tmin_c, temp_c,
     +                          basin_tmax_c, basin_tmin_c,
     +                          basin_tmax_f, basin_tmin_f,
     $                          basin_temp_c, basin_temp_f,
     $                          mru_area,
     +                          basin_area,
     +                          mru_ppt,mru_dep,
     +                          meantmax, meantmin, ppt_lapse, meanppt,
     +                          adjust_snow, tmax_allrain_c,
     +                          tmax_allsnow_c, pptMTH,
     $                          basin_ppt,basin_dep,
     +                          newsnow, mru_rain, mru_snow, pptmix,
     +                          prmx, adjmix_rain,
     +                          tmin_add, tmin_div, tmax_add,
     +                          tmax_div, ppt_add, ppt_div,
     +                          tmin_adj, tmax_adj, rain_code,
     +                          rain_nsta, temp_nuse, rain_nuse,
     +                          temp_meanx, temp_meany, temp_meanz,
     +                          rain_meanx, rain_meany, rain_meanz,
     +                          nrain, tmax_rain_sta, tmin_rain_sta,
     +                          is_rain_day, psta_freq_nuse,
     + solrad_elev, basin_centroid_x, basin_centroid_y, solrad_tmax,
     +                          adjust_rain, nform)

      include 'fmodules.inc'
      
      integer nMRU, temp_nsta, rain_nsta, nform
      integer temp_nuse(MAXTEMP), rain_nuse(MAXRAIN)
      integer psta_freq_nuse(MAXRAIN)
      integer pptmix(MAXMRU), newsnow(MAXMRU)
      integer rain_code(MAXMO), nrain, is_rain_day
      real adjust_snow(MAXMO), adjust_rain(MAXMO)
      real MRUx(MAXMRU), MRUy(MAXMRU), MRUelev(MAXMRU)
      real max_lapse(MAXLAPSE,MAXMO), min_lapse(MAXLAPSE,MAXMO)
      real ppt_lapse(MAXLAPSE,MAXMO)
      real temp_STAx(MAXTEMP),temp_STAy(MAXTEMP), temp_STAelev(MAXTEMP)
      real rain_STAx(MAXRAIN),rain_STAy(MAXRAIN), rain_STAelev(MAXRAIN)
      real tmax_f(MAXMRU), tmin_f(MAXMRU), temp_f(MAXMRU)
      real tmax_c(MAXMRU), tmin_c(MAXMRU), temp_c(MAXMRU)
      real basin_tmax_c, basin_tmin_c, basin_temp_c
      real basin_tmax_f, basin_tmin_f, basin_temp_f
      real mru_area(MAXMRU), basin_area
      real pptMTH(MAXRAIN,MAXMO)
      real meantmax(MAXMO), meantmin(MAXMO)
      real meanppt(MAXMO)
      real mru_ppt(MAXMRU),mru_dep(MAXMRU)
      real tmax_allrain_c(MAXMO), tmax_allsnow_c
      real basin_ppt,basin_dep
      real mru_rain(MAXMRU), mru_snow(MAXMRU), prmx(MAXMRU)
      real adjmix_rain(MAXMO)
      real tmin_add, tmin_div, tmax_add
      real tmax_div, ppt_add, ppt_div
      real tmin_adj(MAXMRU), tmax_adj(MAXMRU)
      real temp_meanx(MAXMO), temp_meany(MAXMO), temp_meanz(MAXMO)
      real rain_meanx(MAXMO), rain_meany(MAXMO), rain_meanz(MAXMO)
      real tmax_rain_sta(MAXRAIN), tmin_rain_sta(MAXRAIN)
      real solrad_elev, basin_centroid_x, basin_centroid_y, solrad_tmax

C Local variables
      integer nowtime(6), iy, im, id
      integer foo
      integer xyz_temp_run, xyz_rain_run

      xyzrun = 1

      call dattim ('now',nowtime)
      id = nowtime(3)
      im = nowtime(2)
      iy = nowtime(1)

      foo = xyz_temp_run (nMRU, MRUx, MRUy, MRUelev,
     +                          max_lapse,
     +                          min_lapse, temp_nsta, temp_STAx,
     +                          temp_STAy,
     +                          temp_STAelev,
     +                          temp_nuse, tmax_f, tmin_f, temp_f,
     +                          tmax_c, tmin_c, temp_c,
     +                          basin_tmax_c, basin_tmin_c,
     +                          basin_tmax_f, basin_tmin_f,
     $                          basin_temp_c, basin_temp_f,
     $                          mru_area,
     +                          basin_area,
     +                          meantmax, meantmin, 
     +                          temp_meanx, temp_meany, temp_meanz,
     +                          tmin_add, tmin_div, tmax_add,
     +                          tmax_div,
     +                          tmin_adj, tmax_adj, im,
     +                          tmax_rain_sta, tmin_rain_sta,
     +                          nrain, rain_STAx, rain_STAy,
     +                          rain_STAelev,
     + solrad_elev, basin_centroid_x, basin_centroid_y, solrad_tmax)

      foo = xyz_rain_run (nMRU, MRUx, MRUy, MRUelev,
     +                          rain_STAx, rain_STAy,
     +                          rain_STAelev,
     +                          rain_nuse, tmax_c, tmin_c,
     +                          mru_area,
     +                          basin_area,
     +                          mru_ppt,mru_dep,
     +                          rain_meanx,
     +                          rain_meany, rain_meanz, ppt_lapse,
     +                          meanppt,
     +                          adjust_snow, tmax_allrain_c,
     +                          tmax_allsnow_c, basin_ppt,basin_dep,
     +                          newsnow, mru_rain, mru_snow, pptmix,
     +                          prmx, adjmix_rain,
     +                          ppt_add, ppt_div,
     +                          rain_code,
     +                          rain_nsta, im, nrain,
     +                          tmax_rain_sta, tmin_rain_sta,
     +                          is_rain_day, psta_freq_nuse,
     +                          adjust_rain, nform)

C      write (*,12)  pptmix
C 12   format(' pptmix ',177i2)
C      write (*,13)  mru_rain
C 13   format(' mru_rain ',300f7.3)

      xyzrun = 0

      return
      end

c***********************************************************************
c
c     Main xyz_dist routine
c

      integer function xyz_dist_web (arg)

      include 'fmodules.inc'

      character(len=*) arg
      CHARACTER(len=256) SVN_ID

      integer xyzdecl, xyzinit, xyzrun, retval


      integer nMRU, temp_nsta, rain_nsta
      integer ntemp, nrain, nform
C      integer nuse(MAXTEMP)
      integer newsnow(MAXMRU), pptmix(MAXMRU)
      integer rain_code(MAXMO), is_rain_day
      integer temp_nuse(MAXTEMP), rain_nuse(MAXRAIN)
      integer psta_freq_nuse(MAXRAIN)
      real adjust_snow(MAXMO), adjust_rain(MAXMO)
      real MRUx(MAXMRU), MRUy(MAXMRU), MRUelev(MAXMRU)
      real max_lapse(MAXLAPSE,MAXMO), min_lapse(MAXLAPSE,MAXMO)
      real ppt_lapse(MAXLAPSE,MAXMO)
      real temp_STAx(MAXTEMP),temp_STAy(MAXTEMP), temp_STAelev(MAXTEMP)
      real rain_STAx(MAXRAIN),rain_STAy(MAXRAIN), rain_STAelev(MAXRAIN)
      real tmax_f(MAXMRU), tmin_f(MAXMRU), temp_f(MAXMRU)
      real tmax_c(MAXMRU), tmin_c(MAXMRU), temp_c(MAXMRU)
      real basin_tmax_c, basin_tmin_c, basin_temp_c
      real basin_tmax_f, basin_tmin_f, basin_temp_f
      real mru_area(MAXMRU), basin_area
      real tmaxMTH (MAXTEMP,MAXMO), tminMTH (MAXTEMP,MAXMO)
      real pptMTH (MAXRAIN,MAXMO)
      real meantmax(MAXMO), meantmin(MAXMO)
      real rain_meanx(MAXMO), rain_meany(MAXMO), rain_meanz(MAXMO)
      real temp_meanx(MAXMO), temp_meany(MAXMO), temp_meanz(MAXMO)
      real meanppt(MAXMO)
      real mru_ppt(MAXMRU),mru_dep(MAXMRU)
      real tmax_allrain_c(MAXMO), tmax_allsnow_c
      real basin_ppt,basin_dep
      real mru_snow(MAXMRU), mru_rain(MAXMRU), prmx(MAXMRU)
      real adjmix_rain(MAXMO)
      real tmin_add, tmin_div, tmax_add, tmax_div
      real ppt_add, ppt_div
      real tmin_adj(MAXMRU), tmax_adj(MAXMRU)
      real tmax_rain_sta(MAXRAIN), tmin_rain_sta(MAXRAIN)
      real solrad_elev, basin_centroid_x, basin_centroid_y, solrad_tmax

      save nMRU, temp_nsta, rain_nsta
      save ntemp, nrain, nform
C      save nuse
      save newsnow, pptmix
      save adjust_snow, adjust_rain
      save MRUx, MRUy, MRUelev
      save max_lapse, min_lapse
      save ppt_lapse
      save temp_STAx, temp_STAy, temp_STAelev
      save rain_STAx, rain_STAy, rain_STAelev
      save tmax_f, tmin_f, temp_f
      save tmax_c, tmin_c, temp_c
      save basin_tmax_c, basin_tmin_c, basin_temp_c
      save basin_tmax_f, basin_tmin_f, basin_temp_f
      save mru_area, basin_area
      save tmaxMTH, tminMTH
      save pptMTH
      save meantmax, meantmin
      save rain_meanx, rain_meany, rain_meanz
      save temp_meanx, temp_meany, temp_meanz
      save meanppt
      save mru_ppt,mru_dep
      save tmax_allrain_c, tmax_allsnow_c
      save basin_ppt,basin_dep
      save mru_snow, mru_rain, prmx
      save adjmix_rain
      save tmin_add, tmin_div, tmax_add, tmax_div
      save ppt_add, ppt_div
      save tmin_adj, tmax_adj
      save rain_code
      save temp_nuse, rain_nuse
      save tmax_rain_sta, tmin_rain_sta
      save is_rain_day
      save psta_freq_nuse
      save solrad_elev, basin_centroid_x, basin_centroid_y, solrad_tmax
      save SVN_ID

      SVN_ID = 
     $     '$Id: xyz_dist_web.f 29 2006-07-06 23:03:45Z rmwebb $ '

      retval = 0

      if(arg.eq.'declare') then
        retval = xyzdecl (tmax_f, tmin_f, temp_f,
     +                           tmax_c, tmin_c, temp_c,
     +                           basin_tmax_c, basin_tmin_c,
     +                           basin_tmax_f, basin_tmin_f,
     $                           basin_temp_c, basin_temp_f,
     +                           mru_ppt, mru_dep,
     $                           basin_ppt, basin_dep,
     $                           newsnow,
     +                           mru_rain, mru_snow, pptmix,
     +                           prmx, tmax_rain_sta, tmin_rain_sta,
     +                           is_rain_day,solrad_tmax)

      else if(arg.eq.'initialize') then

        retval = xyzinit (nMRU, MRUx, MRUy, MRUelev, max_lapse,
     +                           min_lapse, temp_nsta, rain_STAx,
     +                           rain_STAy, rain_STAelev,
     +                           temp_STAx,
     +                           temp_STAy, temp_STAelev,
     +                           mru_area, basin_area,
     +                           tmaxMTH, tminMTH,
     +                           adjust_snow, meantmax, meantmin,
     +                           rain_meanx, rain_meany, rain_meanz,
     +                           temp_meanx, temp_meany, temp_meanz,
     +                           ppt_lapse, meanppt,
     +                           tmax_allrain_c, tmax_allsnow_c,
     +                           pptMTH, adjmix_rain,
     +                           tmin_add, tmin_div, tmax_add,
     +                           tmax_div, ppt_add, ppt_div,
     +                           tmin_adj, tmax_adj,
     +                           rain_code, rain_nsta, ntemp, nrain,
     +                           temp_nuse, rain_nuse, psta_freq_nuse,
     + solrad_elev, basin_centroid_x, basin_centroid_y, adjust_rain,
     + nform)

      else if (arg.eq.'run') then
         retval = xyzrun (nMRU, MRUx, MRUy, MRUelev, max_lapse,
     +                          min_lapse, temp_nsta,
     +                          rain_STAx, rain_STAy, rain_STAelev,
     +                          temp_STAx, temp_STAy, temp_STAelev,
     +                          tmax_f, tmin_f, temp_f,
     +                          tmax_c, tmin_c, temp_c,
     +                          basin_tmax_c, basin_tmin_c,
     +                          basin_tmax_f, basin_tmin_f,
     $                          basin_temp_c, basin_temp_f,
     $                          mru_area,
     +                          basin_area,
     +                          mru_ppt,mru_dep,
     +                          meantmax, meantmin, ppt_lapse, meanppt,
     +                          adjust_snow, tmax_allrain_c,
     +                          tmax_allsnow_c, pptMTH,
     $                          basin_ppt,basin_dep,
     +                          newsnow, mru_rain, mru_snow, pptmix,
     +                          prmx, adjmix_rain,
     +                          tmin_add, tmin_div, tmax_add,
     +                          tmax_div, ppt_add, ppt_div,
     +                          tmin_adj, tmax_adj, rain_code,
     +                          rain_nsta, temp_nuse, rain_nuse,
     +                          temp_meanx, temp_meany, temp_meanz,
     +                          rain_meanx, rain_meany, rain_meanz,
     +                          nrain, tmax_rain_sta, tmin_rain_sta,
     +                          is_rain_day, psta_freq_nuse,
     + solrad_elev, basin_centroid_x, basin_centroid_y, solrad_tmax,
     +                          adjust_rain, nform)

      end if

      xyz_dist_web = retval
      return
      end

c***********************************************************************
c
c     xyz_temp_run - Temperature calculation
c               calculates daily max and min temperature
c               using data from available stations
c               Outputs a daily max and min Temperature by mru elevation
c

C23456789112345678921234567893123456789412345678951234567896123456789712
      integer function xyz_temp_run (nMRU, MRUx, MRUy, MRUelev,
     +                          max_lapse,
     +                          min_lapse, temp_nsta, temp_STAx,
     +                          temp_STAy,
     +                          temp_STAelev,
     +                          temp_nuse, tmax_f, tmin_f, temp_f,
     +                          tmax_c, tmin_c, temp_c,
     +                          basin_tmax_c, basin_tmin_c,
     +                          basin_tmax_f, basin_tmin_f,
     $                          basin_temp_c, basin_temp_f,
     $                          mru_area,
     +                          basin_area,
     +                          meantmax, meantmin, 
     +                          temp_meanx, temp_meany, temp_meanz,
     +                          tmin_add, tmin_div, tmax_add,
     +                          tmax_div,
     +                          tmin_adj, tmax_adj, im,
     +                          tmax_rain_sta, tmin_rain_sta,
     +                          nrain, rain_STAx, rain_STAy,
     +                          rain_STAelev,
     + solrad_elev, basin_centroid_x, basin_centroid_y, solrad_tmax)

      include 'fmodules.inc'
      
      integer nMRU, temp_nsta
      integer temp_nuse(MAXTEMP), im
      integer nrain
      real MRUx(MAXMRU), MRUy(MAXMRU), MRUelev(MAXMRU)
      real max_lapse(MAXLAPSE,MAXMO), min_lapse(MAXLAPSE,MAXMO)
      real temp_STAx(MAXTEMP),temp_STAy(MAXTEMP), temp_STAelev(MAXTEMP)
      real tmax_f(MAXMRU), tmin_f(MAXMRU), temp_f(MAXMRU)
      real tmax_c(MAXMRU), tmin_c(MAXMRU), temp_c(MAXMRU)
      real basin_tmax_c, basin_tmin_c, basin_temp_c
      real basin_tmax_f, basin_tmin_f, basin_temp_f
      real mru_area(MAXMRU), basin_area
      real meantmax(MAXMO), meantmin(MAXMO), temp_meanx(MAXMO)
      real temp_meany(MAXMO), temp_meanz(MAXMO)
      real tmin_add, tmin_div, tmax_add
      real tmax_div
      real tmin_adj(MAXMRU), tmax_adj(MAXMRU)
      real tmax_rain_sta(MAXRAIN), tmin_rain_sta(MAXRAIN)
      real rain_STAx(MAXRAIN), rain_STAy(MAXRAIN), rain_STAelev(MAXRAIN)
      real solrad_elev, basin_centroid_x, basin_centroid_y, solrad_tmax

C Local variables
      integer i, j
      integer ntmin, ntmax
      real tmax(MAXTEMP), tmin(MAXTEMP)
      real intmax, intmin, xmax, xmin, intercept
      real x1,  stmax, stmin
      real sumtmin, sumtmax, xtmin, xtmax, ytmin
      real ytmax, ztmin, ztmax
      real tmax_mru, tmin_mru
      real x_add, x_div, y_add, y_div, z_add, z_div
      real  zrain, xrain, yrain
      real zsolrad, xsolrad, ysolrad

        common /lauren/ x_add, x_div, y_add, y_div, z_add, z_div

      xyz_temp_run = 0

      if(getvar('obs', 'tsta_max_c', MAXTEMP, 'real', tmax)
     +  .ne.0) return 

      if(getvar('xyz_dist', 'tsta_min_c', MAXTEMP, 'real', tmin)
     +  .ne.0) return

      sumtmin=0.0
      sumtmax=0.0
      ntmin=0
      ntmax=0
      xtmin=0.0
      xtmax=0.0
      ytmin=0.0
      ytmax=0.0
      ztmin=0.0
      ztmax=0.0

C Do not Transform the coordinates of the temp stations
c until after summing

      do j=1,temp_nsta
         i = temp_nuse(j)

         if(tmax(i).gt.-55.0)then
            ntmax=ntmax+1
c            sumtmax=sumtmax+((tmax(i) + tmax_add) / tmax_div)
            sumtmax=sumtmax+tmax(i)
            xtmax=xtmax+temp_STAx(i)
            ytmax=ytmax+temp_STAy(i)
            ztmax=ztmax+temp_STAelev(i)
         end if

         if(tmin(i).gt.-55.0)then
            ntmin=ntmin+1
c            sumtmin=sumtmin+((tmin(i) + tmin_add) / tmin_div)
            sumtmin=sumtmin+tmin(i)
            xtmin=xtmin+temp_STAx(i)
            ytmin=ytmin+temp_STAy(i)
            ztmin=ztmin+temp_STAelev(i)
         end if
      end do
c
c calculate means 
c
      if(ntmin.gt.0)then
         stmin=sumtmin/float(ntmin)
         xtmin=xtmin/float(ntmin)
         ytmin=ytmin/float(ntmin)
         ztmin=ztmin/float(ntmin)
        stmin=(stmin + tmin_add) / tmin_div
        xtmin=(xtmin + x_add) / x_div
        ytmin=(ytmin + y_add) / y_div
        ztmin=(ztmin + z_add) / z_div
      else
c these are already transformed
         stmin=meantmin(im)
         xtmin=temp_meanx(im)
         ytmin=temp_meany(im)
         ztmin=temp_meanz(im)
      end if

      if(ntmax.gt.0)then      
         stmax=sumtmax/float(ntmax)
         xtmax=xtmax/float(ntmax)
         ytmax=ytmax/float(ntmax)
         ztmax=ztmax/float(ntmax)
        stmax=(stmax + tmax_add) / tmax_div
        xtmax=(xtmax + x_add) / x_div
        ytmax=(ytmax + y_add) / y_div
        ztmax=(ztmax + z_add) / z_div
      else
c these are already transformed
         stmax=meantmax(im)
         xtmax=temp_meanx(im)
         ytmax=temp_meany(im)
         ztmax=temp_meanz(im)
      end if  

c
c adjust the values if not using all the stations
c
      if(temp_meanz(im).ne.ztmin)then
         intercept=(stmin)-
     .             (min_lapse(3,im)*ztmin)-
     .             (min_lapse(1,im)*xtmin)-
     .             (min_lapse(2,im)*ytmin) 
         stmin=(min_lapse(3,im)*temp_meanz(im))+
     .         (min_lapse(1,im)*temp_meanx(im))+
     .         (min_lapse(2,im)*temp_meany(im))+intercept
      end if
c
      if(temp_meanz(im).ne.ztmax)then
         intercept=(stmax)-
     .             (max_lapse(3,im)*ztmax)-
     .             (max_lapse(1,im)*xtmax)-
     .             (max_lapse(2,im)*ytmax)
         stmax=(max_lapse(3,im)*temp_meanz(im))+
     .         (max_lapse(1,im)*temp_meanx(im))+
     .         (max_lapse(2,im)*temp_meany(im))+intercept
      end if
c

c
c now redistribute based on lapse rates
c                   redistribute to mrus

      xmax=stmax
      xmin=stmin

c      --------------
      intmax=xmax-
     .       (max_lapse(3,im)*temp_meanz(im))-
     .       (max_lapse(1,im)*temp_meanx(im))-
     .       (max_lapse(2,im)*temp_meany(im))
c
      intmin=xmin-
     .       (min_lapse(3,im)*temp_meanz(im))-
     .       (min_lapse(1,im)*temp_meanx(im))-
     .       (min_lapse(2,im)*temp_meany(im))
c

C123456789112345678921234567893123456789412345678951234567896123456789712

      basin_tmax_c = 0.0
      basin_tmin_c = 0.0
      basin_temp_c = 0.0
      basin_tmax_f = 0.0
      basin_tmin_f = 0.0
      basin_temp_f = 0.0

C
C  Compute maximum temperature at XY centroid of basin at the elevation of
C  the solrad station used to develop the DD solrad curves.
C
      zsolrad = (solrad_elev + z_add) / z_div
      xsolrad = (basin_centroid_x + x_add) / x_div
      ysolrad = (basin_centroid_y + y_add) / y_div
      solrad_tmax = (max_lapse(1,im)*xsolrad) +
     .            (max_lapse(2,im)*ysolrad) +
     .            (max_lapse(3,im)*zsolrad) + intmax


      solrad_tmax = (solrad_tmax * tmax_div) - tmax_add

C
C  Compute temperatures at precip stations.
C
      do i=1, nrain
         zrain = (rain_STAelev(i) + z_add) / z_div
         xrain = (rain_STAx(i) + x_add) / x_div
         yrain = (rain_STAy(i) + y_add) / y_div
         tmax_rain_sta(i) =(max_lapse(1,im)*xrain) +
     .            (max_lapse(2,im)*yrain) +
     .            (max_lapse(3,im)*zrain) + intmax

         tmin_rain_sta(i) =(min_lapse(1,im)*xrain)+
     .            (min_lapse(2,im)*yrain)+
     .            (min_lapse(3,im)*zrain)+ intmin

         tmax_rain_sta(i)= (tmax_rain_sta(i) * tmax_div) - tmax_add
         tmin_rain_sta(i)= (tmin_rain_sta(i) * tmin_div) - tmin_add

         if(tmax_rain_sta(i).lt.tmin_rain_sta(i))then
            x1=tmax_rain_sta(i)
            tmax_rain_sta(i)=tmin_rain_sta(i)
            tmin_rain_sta(i)=x1
         end if

      end do

C123456789112345678921234567893123456789412345678951234567896123456789712
      do i=1,nMRU
C         tmax_f(i)=(max_lapse(1,im)*MRUx(i)) +
C     .            (max_lapse(2,im)*MRUy(i)) +
C     .            (max_lapse(3,im)*MRUelev(i)) + intmax
C
C         tmin_f(i)=(min_lapse(1,im)*MRUx(i))+
C     .            (min_lapse(2,im)*MRUy(i))+
C     .            (min_lapse(3,im)*MRUelev(i))+ intmin
C
C
C         tmax_f(i)= (tmax_f(i) * tmax_div) - tmax_add
C         tmin_f(i)= (tmin_f(i) * tmin_div) - tmin_add

C
C  At this point, all temperatures are in the units
C  of the temperatures in the data file.
C
         tmax_mru=(max_lapse(1,im)*MRUx(i)) +
     .            (max_lapse(2,im)*MRUy(i)) +
     .            (max_lapse(3,im)*MRUelev(i)) + intmax

         tmin_mru=(min_lapse(1,im)*MRUx(i))+
     .            (min_lapse(2,im)*MRUy(i))+
     .            (min_lapse(3,im)*MRUelev(i))+ intmin

C
C  Transform back
C
         tmax_mru = (tmax_mru * tmax_div) - tmax_add
         tmin_mru = (tmin_mru * tmin_div) - tmin_add

C
C  Temp adjustment by MRU
C
         tmax_mru = tmax_mru + tmax_adj(i)
         tmin_mru = tmin_mru + tmin_adj(i)

C
C  If max is less than min, switch
C
         if (tmax_mru .lt. tmin_mru) then
            x1 = tmax_mru
            tmax_mru = tmin_mru
            tmin_mru = x1
         end if

C
C  Now sort out units.
c
c  Not needed since temperature inputs are always celsius - RW
C

            tmax_c(i) = tmax_mru
            tmin_c(i) = tmin_mru
            temp_c(i) = (tmax_c(i) + tmin_c(i)) / 2.0
            tmax_f(i) = (tmax_c(i) * (9.0 / 5.0)) + 32.0
            tmin_f(i) = (tmin_c(i) * (9.0 / 5.0)) + 32.0
            temp_f(i) = (tmax_f(i) + tmin_f(i)) / 2.0

            basin_tmax_c = basin_tmax_c + (tmax_c(i) * mru_area(i))
            basin_tmin_c = basin_tmin_c + (tmin_c(i) * mru_area(i))
            basin_temp_c = basin_temp_c + (temp_c(i) * mru_area(i))
   
         end do

      basin_tmax_c = basin_tmax_c / basin_area
      basin_tmin_c = basin_tmin_c / basin_area
      basin_temp_c = basin_temp_c / basin_area
      basin_tmax_f = basin_tmax_c*1.8+32.
      basin_tmin_f = basin_tmin_c*1.8+32.
      basin_temp_f = basin_temp_c*1.8+32.

      return
      end


C23456789112345678921234567893123456789412345678951234567896123456789712
      integer function xyz_rain_run (nMRU, MRUx, MRUy, MRUelev,
     +                          rain_STAx, rain_STAy,
     +                          rain_STAelev,
     +                          rain_nuse, tmax_c, tmin_c,
     +                          mru_area,
     +                          basin_area,
     +                          mru_ppt,mru_dep,
     +                          rain_meanx,
     +                          rain_meany, rain_meanz, ppt_lapse,
     +                          meanppt,
     +                          adjust_snow, tmax_allrain_c,
     +                          tmax_allsnow_c,
     $                          basin_ppt,basin_dep,
     +                          newsnow, mru_rain, mru_snow, pptmix,
     +                          prmx, adjmix_rain,
     +                          ppt_add, ppt_div,
     +                          rain_code,
     +                          rain_nsta, im, nrain,
     +                          tmax_rain_sta, tmin_rain_sta,
     +                          is_rain_day, psta_freq_nuse,
     +                          adjust_rain, nform)

      include 'fmodules.inc'
      
      integer nMRU, rain_nsta, nform

      integer rain_nuse(MAXRAIN)
      integer psta_freq_nuse(MAXRAIN)
      integer newsnow(MAXMRU), pptmix(MAXMRU)
      integer rain_code(MAXMO)
      integer im, nrain, is_rain_day
      integer rain_day
      real adjust_snow(MAXMO), adjust_rain(MAXMO), adjmix
      real MRUx(MAXMRU), MRUy(MAXMRU), MRUelev(MAXMRU)
      real ppt_lapse(MAXLAPSE,MAXMO)
      real rain_STAx(MAXRAIN),rain_STAy(MAXRAIN)
      real rain_STAelev(MAXRAIN)
      real tmax_c(MAXMRU), tmin_c(MAXMRU)
      real mru_area(MAXMRU), basin_area
      real rain_meanx(MAXMO)
      real rain_meany(MAXMO), rain_meanz(MAXMO)
      real meanppt(MAXMO)
      real mru_ppt(MAXMRU),mru_dep(maxmru)
      real tmax_allrain_c(MAXMO), tmax_allsnow_c
      real basin_ppt,basin_dep
      real mru_rain(MAXMRU), mru_snow(MAXMRU), prmx(MAXMRU)
      real adjmix_rain(MAXMO)
      real ppt_add, ppt_div
      real tmax_rain_sta(MAXRAIN)
      real tmin_rain_sta(MAXRAIN)


C Local variables
      integer i, j, iflag, err_chk
      integer nppt
      integer nsta_used
      integer form_data(MAXFORM)
      real ppt(MAXRAIN)

      real intppt, intercept
      real sppt
      real sumppt, xppt
      real yppt, zppt
      real xsumppt
      real x_add, x_div, y_add, y_div, z_add, z_div

        common /lauren/ x_add, x_div, y_add, y_div, z_add, z_div

C DANGER
      integer nowtime(6), iy, id

      xyz_rain_run = 0

      if(getvar('obs', 'precip', MAXRAIN, 'real', ppt)
     +  .ne.0) return

C23456789112345678921234567893123456789412345678951234567896123456789712

c      if(get var('xyz_dist', 'form_data', MAXFORM, 'integer', form_data)
c     +   .ne.0) return

c      form_data(1) = 0

      if (nform .ne. 0) then
         if(getvar('obs', 'form_data', MAXFORM, 'integer',
     +      form_data)
     +      .ne.0) return
      else
         form_data(1) = 0
      endif

c      write (*,*) id, form_data(1)

      if(getvar('xyz_dist', 'rain_day', 1, 'integer', rain_day)
     +   .ne.0) return

C
C Code to check the rain_code parameter to determine if it is
C raining in the basin.
C
      is_rain_day=0
      if (rain_code(im) .eq. 1) then
         do j=1, rain_nsta
            i = rain_nuse(j)
            if(ppt(i).gt.0.0)is_rain_day=1
         end do

      else if (rain_code(im) .eq. 2) then
         do i=1,nrain
            if(ppt(i).gt.0.0)is_rain_day=1
         end do

      else if (rain_code(im) .eq. 3) then
         is_rain_day=1

      else if (rain_code(im) .eq. 4) then
         if (rain_day .eq. 1) then
            is_rain_day=1
         end if

      else if (rain_code(im) .eq. 5) then
         do i=1, nrain
            if(psta_freq_nuse(i).eq.1) then
               if(ppt(i).gt.0.0)is_rain_day=1
            end if
         end do

      end if

C DANGER write
C      call dattim ('now',nowtime)
C      id = nowtime(3)
C      iy = nowtime(1)
C      write(*,*)id,im,iy, '   is_rain_day = ', is_rain_day,
C     +           '     uncorrected CB ppt = ',ppt(6)

c
c add adjust_snow and adjust_rain here
c
      if (rain_code(im) .eq. 1) then
         nsta_used = rain_nsta
      else
         nsta_used = nrain
      end if

      do j = 1, nsta_used
         if (rain_code(im) .eq. 1) then
            i = rain_nuse(j)
         else
            i = j
         end if

         if(ppt(i).lt.0.0)go to 333

         err_chk=3     ! Assume mixed rain snow
         if(tmax_rain_sta(i).le.tmax_allsnow_c)then
            err_chk=1
         else if(tmin_rain_sta(i).gt.tmax_allsnow_c.or.
     .          tmax_rain_sta(i).ge.tmax_allrain_c(im))then
            err_chk=0
         else                ! mixed so use same logic as in precip_web
c                              to weight the adjustment factors.
            adjmix =  ((tmax_rain_sta(i)-tmax_allsnow_c)/
     $           (tmax_rain_sta(i)-tmin_rain_sta(i)))
            if(adjmix.gt.1) adjmix = 1.0
            adjmix = adjmix * adjust_rain(im) +
     $           (1.0-adjmix)*adjust_snow(im)
         end if

         if (err_chk.eq.1) then
            ppt(i)=(ppt(i)*adjust_snow(im))+ppt(i)
         else if (err_chk.eq.0) then
            ppt(i)=(ppt(i)*adjust_rain(im))+ppt(i)
         else
            ppt(i)=(ppt(i)*adjmix)+ppt(i)
         end if

C DANGER write
C         if (j .eq. 6) then
C            write(*,*) "Crested Butte corrected precip = ", ppt(i)
C         end if

333      continue
      end do

C      write(*,*)id,im,iy, ppt
C      write(*,*) '     corrected CB precip = ', ppt(6)

      sumppt=0.0
      nppt=0
      xppt=0.0
      yppt=0.0
      zppt=0.0

      do j=1,rain_nsta
         i = rain_nuse(j)
         if(ppt(i).ge.0.0)then
            nppt=nppt+1
c            sumppt=sumppt+((ppt(i) + ppt_add) / ppt_div)
            sumppt=sumppt+ppt(i)
            xppt=xppt+rain_STAx(i)
            yppt=yppt+rain_STAy(i)
            zppt=zppt+rain_STAelev(i)
         end if
      end do
c
c calculate means 
c
      if(nppt.gt.0)then      
         sppt=sumppt/float(nppt)
         xppt=xppt/float(nppt)
         yppt=yppt/float(nppt)
         zppt=zppt/float(nppt)
        sppt=(sppt + ppt_add) / ppt_div
        xppt=(xppt + x_add) / x_div
        yppt=(yppt + y_add) / y_div
        zppt=(zppt + z_add) / z_div
      else
c these are already transformed
         sppt=meanppt(im)
         xppt=rain_meanx(im)
         yppt=rain_meany(im)
         zppt=rain_meanz(im)
      end if  

C      call dattim ('now',nowtime)
C      id = nowtime(3)
C      iy = nowtime(1)
C      write(*,*)id,im,iy, sppt, xppt, yppt, zppt
c
c adjust the values if not using all the stations
c
      if(rain_meanz(im).ne.zppt)then
c         if(iflag.eq.0)then
C            print *,'iflag=',iflag
            intercept=(sppt)-
     .                (ppt_lapse(3,im)*zppt)-
     .                (ppt_lapse(1,im)*xppt)-
     .                (ppt_lapse(2,im)*yppt)
            sppt=(ppt_lapse(3,im)*rain_meanz(im))+
     .           (ppt_lapse(1,im)*rain_meanx(im))+
     .           (ppt_lapse(2,im)*rain_meany(im))+intercept

c        end if
      end if

C      write(*,*)id,im,iy, intercept, sppt

      xppt= sppt

      intppt=xppt-
     .       (ppt_lapse(3,im)*rain_meanz(im))-
     .       (ppt_lapse(1,im)*rain_meanx(im))-
     .       (ppt_lapse(2,im)*rain_meany(im))

      basin_ppt = 0.0
      basin_dep = 0.0
C      write(*,*)       ' rain_meanz(im) = ', rain_meanz(im),
C     + '   intppt = ', intppt

C123456789112345678921234567893123456789412345678951234567896123456789712
      do i=1,nMRU
         newsnow(i) = 0

         if(is_rain_day.eq.0)then
            mru_ppt(i)=0.0
         else
            mru_ppt(i)=(ppt_lapse(1,im)*MRUx(i))+
     .                 (ppt_lapse(2,im)*MRUy(i))+
     .                 (ppt_lapse(3,im)*MRUelev(i)) + intppt

            mru_ppt(i)= mru_ppt(i) * ppt_div - ppt_add

            if(mru_ppt(i).lt.0.0)mru_ppt(i)=0.0
         end if

C******Zero precipitation on MRU

         if (mru_ppt(i) .eq. 0) then
            pptmix(i) = 0
            mru_rain(i) = 0.
            mru_snow(i) = 0.
            prmx(i) = -0.1
            goto 1000
         end if


C******If observed temperature data are not available or if observed
C******form data are available and rain is explicitly specified then
C******precipitation is all rain.

         if(form_data(1).eq.2) then
            pptmix(i) = 0
            mru_rain(i) = mru_ppt(i)
            mru_snow(i) = 0.
            prmx(i) = 1.0

C******If form data are available and snow is explicitly specified or if
C******maximum temperature is below the base temperature for snow then
C******precipitation is all snow

         else if(form_data(1).eq.1.or.tmax_c(i).le.tmax_allsnow_c) then
            pptmix(i) = 0
            mru_rain(i) = 0.
            mru_snow(i) = mru_ppt(i)
            newsnow(i) = 1
            prmx(i) = 0.0

C******If minimum temperature is above base temperature for snow or
C******maximum temperature is above all_rain temperature then
C******precipitation is all rain

         else if(tmin_c(i).ge.tmax_allsnow_c.or.
     +          tmax_c(i).ge.tmax_allrain_c(im)) then
            pptmix(i) = 0
            mru_rain(i) = mru_ppt(i)
            mru_snow(i) = 0.
            prmx(i) = 1.0

C******Otherwise precipitation is a mixture of rain and snow

         else
            prmx(i) = ((tmax_c(i)-tmax_allsnow_c)/
     $           (tmax_c(i)-tmin_c(i)))* adjmix_rain(im)

C******Unless mixture adjustment raises the proportion of rain to
C******greater than or equal to 1.0 in which case it all rain

            if(prmx(i).ge.1.) then
               pptmix(i) = 0
               prmx(i) = 1.0
               mru_rain(i) = mru_ppt(i)
               mru_snow(i) = 0.

C******If not, it is a rain/snow mixture

            else
c
c           RMTW - the following logic is OK here (changed in precip_web)
c                  since the initial mixed snow/rain correction was applied
c                  above.
c
c            mru_ppt(i) = ppt * pcor
c            mru_rain(i) = prmx(i) * mru_ppt(i)
c            mru_snow(i) = mru_ppt(i) - mru_rain(i)
               pptmix(i) = 1
               mru_rain(i) = prmx(i) * mru_ppt(i)
               mru_snow(i) = mru_ppt(i) - mru_rain(i)
               newsnow(i) = 1
            end if
         end if

 1000    continue
         mru_dep(i) = mru_ppt(i)
         basin_ppt = basin_ppt + mru_ppt(i) * mru_area(i)
      end do

      basin_ppt = basin_ppt/basin_area
      basin_dep = basin_ppt
      return
      end

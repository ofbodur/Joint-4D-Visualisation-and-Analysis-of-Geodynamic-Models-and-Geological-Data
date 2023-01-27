#!/usr/bin/python
#=====================================================================
#                      Python Scripts for CitComS 
#         Preprocessing, Data Assimilation, and Postprocessing
#
#                    AUTHORS: Mark Turner, Dan Bower
#
#                  ---------------------------------
#             (c) California Institute of Technology 2007
#                        ALL RIGHTS RESERVED
#=====================================================================
'''A set of general purpose functions for working with Citcom Data.

This module has functions to read and parse various file types associated 
with CitComS and closely related data sets.

Most functions take a dictionary of arguments to pass and adjust paramters.
'''
#=====================================================================
#=====================================================================
# imports
# import getopt, os, string, sys, math, time, commands # python2 version
import getopt, os, string, sys, math, time, subprocess # python3 version, replaced commands with subprocess
import datetime
import pprint

import CitcomParser, Core_Util

from Core_Util import now
#=====================================================================
#=====================================================================
# global variables
verbose = False

r2d = 180.0 / math.pi
earth_radius = 6371.0

#=====================================================================
#=====================================================================
def parse_citcoms_cfg_file( filename ):
    '''Parse a 'key = value' style control file, and return a dictionary of those settings;

Required arguments:
 filename of pid file to parse
    
Optional arguments:
 none

Output arguments:
 none

Return value:
 none
'''

    # empty settings dict
    settings = {}

    file = open(filename)
    if verbose: print (now(), 'parse_citcoms_cfg_file: read: %(filename)s') % vars()

    try : 
        # read into a list of lines
        lines = file.read().splitlines()
    finally:
        file.close()

    # read
    for line in lines:

        # skip blank lines
        # if string.strip(line) == '': # python2 version.
        # print(line.strip())
        if line.strip() == '': # python3 version, replaced strip with line.strip()
            continue # to next line in control file

        # skip lines starting with '#'
        if line.startswith('#'):
            continue # to next line in control file

        # process section boundary into new key 
        if line.startswith('[') :
            continue # to next line in control file

        # key/value pairs can be seperated by whitespaces
        opt = line.rstrip()
        opt = opt.replace(' ','') # remove any spaces padding

        # isolate the key and value
        # (key,val) = string.split(opt, '=') # python2 version
        # print(opt)
        # print(opt.split('='),'This Line')
        key = opt.split('=')[0] # python3 version
        # print(key, "is the key")
        val = opt.split('=')[1]  # python3 version
        # WARN: this can be really verbose:
        # if verbose:  
        #   print now(), 'parse_citcoms_cfg_file: key , val = ',key,',',val

        # this key val pair is a general setting
        settings[key] = val

    # WARN: this can be really verbose:
    #if verbose:
    #    print now(), 'parse_citcoms_cfg_file: settings ='
    #    pprint.PrettyPrinter(indent=2).pprint(settings)

    return settings
#=====================================================================
#=====================================================================
def set_region( args ):
    '''from citcoms parameters, determine region

Required arguments:
 args['citcoms_pid'] = name of pid file of model
 args['nproc_surf'] = used to determin model type 

Optional arguments:
 none

Output arguments:
 args['zone'] = 'regional' , or 'global'
 args['latmax'] = bounds of model 
 args['latmin'] = bounds of model 
 args['lonmin'] = bounds of model 
 args['lonmax'] = bounds of model 
 args['R'] = bounds of model 

Return value:
 none
''' 


    # parse the parameter file; CitcomParser checks for existance
    d = parse_citcoms_cfg_file( args['citcoms_pid'] )
    args.update(d)

    # determine type of CitComS run:
    nproc_surf = int( args['nproc_surf'] )

    if nproc_surf == 1:
        args['zone'] = 'regional'

        args['latmax'] = 90.0 - ( r2d * float(args['theta_min']) )
        args['latmin'] = 90.0 - ( r2d * float(args['theta_max']) )

        args['lonmin'] = ( r2d * float(args['fi_min']) )
        args['lonmax'] = ( r2d * float(args['fi_max']) )

        # correct bounds 
        if args.get('force360'):
            if args['lonmin'] < 0:
                args['lonmin'] = args['lonmin'] + 360
            if args['lonmax'] < 0:
                args['lonmax'] = args['lonmax'] + 360

        # correct bounds 
        if args.get('force180'):
            if args['lonmin'] > 180:
                args['lonmin'] = args['lonmin'] - 360
            if args['lonmax'] > 180:
                args['lonmax'] = args['lonmax'] - 360

        args['R'] = "%g/%g/%g/%g" % ( \
            args['lonmin'], args['lonmax'], \
            args['latmin'], args['latmax'] )

    elif nproc_surf == 12:
        args['zone'] = 'global'
        args['R'] = '0/360/-90/90'

    else:
        args['zone'] = None
        msg = 'Unable to determine model type: gloabal or regional.'
        raise (ValueError, msg)

    # re-set the zone to 2D for 'thin' models 
    if args['nodex'] <= 3 or args['nodey'] <= 3 :
        args['zone'] = '2D'
        args['lonmin'] -= 1.0
        args['lonmax'] += 1.0
        args['R'] = "%g/%g/%g/%g" % ( \
            args['lonmin'], args['lonmax'], \
            args['latmin'], args['latmax'] )


    if verbose: 
        print (now(), 'set_model_region: datafile =', args['datafile'])
        print (now(), 'set_model_region: zone =', args['zone'])
        print (now(), 'set_model_region: R =', args['R'])
#=====================================================================
#=====================================================================
#=====================================================================
def write_cap_polygons(args):
    '''Read a set of Citcoms cap files, write a single .xy files
holding the outline of the cap boundaries.

Required arguments:
 args['citcoms_pid'] = name of pid file;
 args['time'] = output time in model steps;

Optional arguments:
 args['cap_list'] = list of cap numbers [0, 11] as eithern ints or strings
  default is to read all caps: 0 to 11 for global model; 0 for regional model.

Output arguments:
 none

Return value:
 none

WARNING: 
WARNING: work in progres, does not work for all caps yet
WARNING: 

'''

    # parse citcom model data
    citcoms = parse_citcoms_cfg_file( args['citcoms_pid'] )
    datafile = citcoms['datafile']
    nodex = int( citcoms['nodex'] )
    nodey = int( citcoms['nodey'] )
    nodez = int( citcoms['nodez'] )
    nproc_surf = int( citcoms['nproc_surf'] )

    time = args['time']

    # FIX: do we need level?
    z = 0 

    # set cap list from nproc_surf
    cap_list = [0]
    if nproc_surf == 12:
        cap_list = range(12)

    # reset cap list if requested 
    if args.get('cap_list'):
        cap_list = args.get('cap_list')

    if verbose: print (now(), 'write_cap_polygons: cap_list ', cap_list)

    # final output file name
    poly_file = 'polygon.cap_bounds.xy' % vars()
    out_poly= open(poly_file, 'w')
    if verbose: print (now(), 'write_cap_polygons: open(w) ', poly_file)

    gmtchar=">"

    # loop over caps
    for c in cap_list:
        c = '%02d' % int(c)

        # open the cap input file
        inp = '%(datafile)s.cap%(c)s.%(time)s' % vars()
        input = open(inp)

        nodex, nodey, nodez = input.readline().split('x')
        nodex = int(nodex)
        nodey = int(nodey)
        nodez = int(nodez)

        lat_deg=range(nodex*nodey)
        lon_deg=range(nodex*nodey)

        if verbose: 
            print (now(), 'write_cap_polygons: loop over file: ', inp)
            print (now(), 'write_cap_polygons: loop over j: 0,', nodey)
            print (now(), 'write_cap_polygons: loop over i: 0,', nodex)
            print (now(), 'write_cap_polygons: loop over k: 0,', nodez)

        for i in range(nodey):
            for j in range(nodex): 
                for k in range(nodez):
                    line = input.readline()
                    n=i*nodex+j
                    if k == z :
                        lat = float(line.split(' ')[0])
                        lon = float(line.split(' ')[1])
                        lat_deg[n] = 90.0 - lat * r2d
                        lon_deg[n] = lon * r2d
                        #if int(c) >= 8:
                        #    if lon_deg[n] > 180.0:
                        #        lon_deg[n] = lon_deg[n]-360.0

        input.close()

        #March around the 4 sides of the polygon
        #out_poly.write("%s\n" % gmtchar )
        for i in range(nodex):
            j=0 
            n=i*nodey+j
            out_poly.write("%g %g\n" % (lon_deg[n],lat_deg[n]) )
    
        #out_poly.write("%s\n" % gmtchar )
        for j in range(nodey):
            i=nodex-1 
            n=i*nodey+j
            out_poly.write("%g %g\n" % (lon_deg[n],lat_deg[n]) )

        #out_poly.write("%s\n" % gmtchar )
        for i in range(nodex):
            j=nodey-1
            n=(nodex-1-i)*nodey+j
            out_poly.write("%g %g\n" % (lon_deg[n],lat_deg[n]) )

        #out_poly.write("%s\n" % gmtchar )
        for j in range(nodey):
            i=0 
            n=i*nodey+(nodey-1-j)
            out_poly.write("%g %g\n" % (lon_deg[n],lat_deg[n]) )

        # separate each cap boundary
        out_poly.write('>\n')

    out_poly.close()

    return poly_file
#=====================================================================
#=====================================================================
def slice( args ):
    '''Create a 2-D slice from a set of Citcoms cap files.
    Slice type may be one of:
    'zslice' for fixed radius, with data written in 'lat lon data',
    'yslice' for fixed lon, with data written in 'lat radius data',
    'xslice' for fixed lat, with data written in 'lon radius data',
    where lat ranges -90 to 90 degrees, lon ranges 0 to 360 degrees, radius ranges 1.0 to 0.55 non-dimensional values;

Required arguments:
 args['slice_type'] = type of slice to perform

Optional arguments:
 none

Output arguments:
 none

Return value:
 none
'''

    type = args.get('slice_type')

    if   type == 'zslice': zslice( args )
    elif type == 'yslice': yslice( args ) 
    elif type == 'xslice': xslice( args ) 
    else:
        msg = '''Unknown slice type;
Slice type may be one of:
    'zslice' for a fixed radius, with data written in 'lat lon data',
    'yslice' for a fixed lon, with data written in 'lat radius data',
    'xslice' for a fixed lat, with data written in 'lon radius data',
where lat ranges -90 to 90, lon ranges 0 to 360, radius ranges 1.0 to 0.55'''
        raise (IndexError, msg)
#=====================================================================
#=====================================================================
def zslice( args = {} ):
    '''Read a set of Citcoms cap files, write a set of .xyz files,
and return the set of file names as a comma separated list.
Output file name format: datafile.time.z.field_id
 One output .xyz file will be produced for each field_id given.

Required arguments:
 args['citcoms_pid'] = name of pid file;
 args['time'] = model time steps [0,model end step];
 args['z'] = model node number in radial direction: [0, nodez - 1];
 args['field'] = comma separated list of field ids: 
 {temp, visc, velocity} (NOTE velocity parallel to slice plane);

Optional arguments: 
 args['velocity_scale'] = float value to scale velocity magnitude
NOTE: velocity data is produced in cm/yr scaled by this factor.
 args['velocity_increment'] = float value, <= 1.0, to decimate data

 args['cap_list'] = used to explicitly set which caps to read and write

 args.get('get_lat_min') = set limits on output region
 args.get('get_lat_max') = set limits on output region
 args.get('get_lon_min') = set limits on output region
 args.get('get_lon_max') = set limits on output region

Output arguments:
 none

Return value:
 list of output files written
'''

    # parse citcom model data
    citcoms = parse_citcoms_cfg_file( args['citcoms_pid'] ) 
    datafile = citcoms['datafile']
    nodex = int( citcoms['nodex'] ) 
    nodey = int( citcoms['nodey'] ) 
    nodez = int( citcoms['nodez'] ) 
    nproc_surf = int( citcoms['nproc_surf'] ) 

    # set cap list from nproc_surf
    cap_list = [0] 
    if nproc_surf == 12: 
        cap_list = range(12)

    # reset cap list if requested 
    if args.get('cap_list'):
        cap_list = args.get('cap_list')

    # parse time value
    time = args['time']

    # parse z value
    z = args['z']

    # parse field value, and settings
    value = args['field'] 
    # create a list with the value
    field_list = [ value ]
    # check for comma separated multiple values
    if value.count(','):
        # update list
        field_list = value.split(',')
    if verbose: print (now(), 'zslice: field_list ', field_list)

    # parse velocity scale
    velocity_scale = 1.0
    if args.get('velocity_scale'):
        velocity_scale = float(args.get('velocity_scale'))

    # scale velocities in cm/yr 
    thermdiff = float( citcoms['thermdiff'] )
    layer_km = earth_radius
    base_scale = (thermdiff/(layer_km*1e3))*(100.0*3600.0*24.0*365.25)

    # parse velocity increment 
    velocity_increment = 1.0
    if args.get('velocity_increment'): 
        velocity_increment = float(args.get('velocity_increment'))

    # parse geo bounds
    check_geographic_bounds = False
    if args.get('get_lat_min') \
    or args.get('get_lat_max') \
    or args.get('get_lon_min') \
    or args.get('get_lon_max') :
        check_geographic_bounds = True
 
    # empty lists to hold data lines
    vx_list = []
    vy_list = []
    vz_list = []
    temp_list = []
    visc_list = []
    composition_list = []
    velocity_list = []

    # loop over caps
    for c in cap_list:
        c = '%02d' % int(c)

        # open the cap input file 
        inp = '%(datafile)s.cap%(c)s.%(time)s' % vars()
        inp_file = open(inp)
        if verbose: print (now(), 'zslice: open ', inp)

        # read the cap file into one big list of lines
        lines = inp_file.readlines()
        inp_file.close()
        if verbose: print (now(), 'zslice: close ', inp)

        # remove the header line
        lines = lines[1:]

        # composition file
        if 'comp' in field_list:
            print ("opening composition file")
            # form file name from model parameters
            opt = '%(datafile)s.opt%(c)s.%(time)s' % vars()
            # open the file 
            opt_file = open(opt)
            opt_lines = opt_file.readlines()
            # pop the header line
            opt_lines = opt_lines[1:]

        if verbose: print (now(), 'zslice: loop over file: ', inp)
        if verbose: print (now(), 'zslice: loop over j: 0,', nodey)
        if verbose: print (now(), 'zslice: loop over i: 0,', nodex)
        if verbose: print (now(), 'zslice: loop over k: 0,', nodez)
        # loop over the node numbers
        # NOTE loop order: Y , X , Z 
        for j in range(nodey):
            for i in range(nodex):
                for k in range(nodez):

                    # process only that data at level z
                    if not k == int(args['z']):
                        continue # to next k

                    # compute line number
                    n = k + i*nodez + j*nodez*nodex

                    line = lines[n].rstrip()

                    # separate into cols: lat lon vx vy vz temp visc 
                    cols = line.split(' ')

                    # convert angles
                    lat = 90 - float(cols[0]) * r2d 
                    lon =      float(cols[1]) * r2d 
                    # non-dimensional radius
                    rad =      float(cols[2]) 
                    vx  =      float(cols[3])
                    vy  =      float(cols[4])
                    vz  =      float(cols[5])
                    temp =     float(cols[6])
                    visc =     math.log10(float(cols[7])) # log scale

                    if 'comp' in field_list:
                        opt_line = opt_lines[n].rstrip()

                    if check_geographic_bounds:
                        if lat < float(args['get_lat_min']) \
                        or lat > float(args['get_lat_max']) \
                        or lon < float(args['get_lon_min']) \
                        or lon > float(args['get_lon_max']) :
                            continue # to next k

                    # assemble data lists
                    temp_list.append( '%f %f %f' % (lon, lat, temp) )
                    visc_list.append( '%f %f %f' % (lon, lat, visc) )

                    if 'comp' in field_list:
                        composition_list.append( '%f %f %s' % (lon, lat,opt_line) )

                    # check for combined fields 
                    if field_list.count('velocity'):
                        if vy <= 0.0:
                              if vx < 0.0:
                                  a = math.degrees(math.atan(vy/vx))
                                  azimuth = 360.0 - a
                              elif vx > 0.0:
                                  a = math.degrees(math.atan(-vy/vx))
                                  azimuth = 180.0 + a
                              else:
                                  azimuth = 270.0
                        elif vy > 0.0:
                              if vx < 0.0:
                                  a = math.degrees(math.atan(-vy/vx))
                                  azimuth = a
                              elif vx > 0.0:
                                  a = math.degrees(math.atan(vy/vx))
                                  azimuth = 180.0 - a
                              else:
                                  azimuth = 90.0
    
                        # apply scale 
                        length = (base_scale * math.hypot(vx,vy)) / velocity_scale
                        data = '%f %f %f %f' % (lon, lat, azimuth, length)
                        velocity_list.append( data ) 
                    # end of velocity calcs

                # end of k loop over Z
            # end of i loop over X
        # end of j loop over Y
    # end of loop over caps



    # write out data files
    out_list = [] 

    out_prefix = datafile + '.' \
               + args['time'] + '.' \
               + args['z'] + '.'

    if field_list.count('temp'): 
        # open output file
        out_name = out_prefix + 'temp' + '.xyz'
        out_file = open(out_name, 'w')
        if verbose: print (now(), 'zslice: open,w ', out_name)
        for line in temp_list:
            out_file.write( '%s\n' % line )
        out_file.close()
        out_list.append(out_name)
        if verbose: print (now(), 'zslice: close: ', out_name)

    if field_list.count('visc'): 
        # open output file
        out_name = out_prefix + 'visc' + '.xyz'
        out_file = open(out_name, 'w')
        if verbose: print (now(), 'zslice: open,w ', out_name)
        for line in visc_list:
            out_file.write( '%s\n' % line )
        out_file.close()
        out_list.append(out_name)
        if verbose: print (now(), 'zslice: close: ', out_name)

    if 'comp' in field_list: 
        out_name = out_prefix + 'comp' + '.xyz'
        out_file = open(out_name, 'w')
        if verbose: print (now(), 'zslice: open,w ', out_name)
        for line in composition_list:
            out_file.write( '%s\n' % line )
        out_file.close()
        if verbose: print (now(), 'zslice: close: ', out_name)

    if field_list.count('velocity'): 
        # open output file
        out_name = out_prefix + 'velocity' + '.xyz'
        out_file = open(out_name, 'w')
        if verbose: print (now(), 'zslice: open,w ', out_name)
        p = 0
        q = 0 
        for line in velocity_list:
            # check for sub-sample
            q = q + 1
            if q % int(1/float(velocity_increment)) != 0:
                continue # to next line
            out_file.write( '%s\n' % line )
        out_file.close()
        out_list.append(out_name)
        if verbose: print (now(), 'zslice: close: ', out_name)

    return out_list
#=====================================================================
#=====================================================================
#=====================================================================
def xslice(args):
    '''Read a set of Citcoms cap files, write a set of .xyz files,
and return the set of file names as a comma separated list.
Output file name format: datafile.time.x.field_id
One output .xyz file will be produced for each field_id given.


Required arguments:
 args['citcoms_pid'] = name of pid file;
 args['time'] = model time steps [0,model end step];
 args['x'] = model node number in theta direction: [0, nodex - 1];
 args['field'] = comma separated list of field ids: 
 {temp, visc, velocity, comp} (NOTE velocity parallel to slice plane);

Optional arguments: 
 args['velocity_scale'] = float value to scale velocity magnitude
NOTE: velocity data is produced in cm/yr multiplied by this factor.
 args['velocity_increment'] = float value, <= 1.0, to decimate data
 
 args['cap_list'] = used to explicitly set which caps to read and write

 args.get('get_r_min') = set limits on output region
 args.get('get_r_max') = set limits on output region
 args.get('get_lon_min') = set limits on output region
 args.get('get_lon_max') = set limits on output region

 args['polar_angular_offset'] = angle in degrees 
  Setting this value will cause the velocity vectors' azimuth value
  to be adjusted such that an azimuth of 0 will point 'locally up'.
  This is used to match the polar projection with GMT's PSXY.

Output arguments:
 none

Return value:
 list of output files written 
'''

    # parse citcom model data
    citcoms = parse_citcoms_cfg_file( args['citcoms_pid'] ) 
    datafile = citcoms['datafile']
    nodex = int( citcoms['nodex'] ) 
    nodey = int( citcoms['nodey'] ) 
    nodez = int( citcoms['nodez'] ) 
    nproc_surf = int( citcoms['nproc_surf'] ) 

    # set cap list from nproc_surf
    cap_list = [0] 
    if nproc_surf == 12: 
        cap_list = range(12)

    # reset cap list if requested 
    if args.get('cap_list'):
        cap_list = args.get('cap_list')

    # parse time value
    time = args['time']

    # parse x value
    x = args['x']

    # parse field value, and settings
    value = args['field'] 
    # create a list with the value
    field_list = [ value ]
    # check for comma separated multiple values
    if value.count(','):
        # update list
        field_list = value.split(',')
    if verbose: print (now(), 'xslice: field_list ', field_list)

    # parse velocity scale
    velocity_scale = 1.0
    if args.get('velocity_scale'):
        velocity_scale = float(args.get('velocity_scale'))

    # convert velocities to cm/yr
    thermdiff = float( citcoms['thermdiff'] )
    layer_km = earth_radius
    base_scale = (thermdiff/(layer_km*1e3))*(100.0*3600.0*24.0*365.25)

    # parse velocity increment 
    velocity_increment = 1.0
    if args.get('velocity_increment'): 
        velocity_increment = float(args.get('velocity_increment'))
  
    # phase controls
    phase_args = args
    phase_args.update( citcoms )

    # parse geo bounds
    check_geographic_bounds = False
    if args.get('get_r_min') \
    or args.get('get_r_max') \
    or args.get('get_lon_min') \
    or args.get('get_lon_max') :
        check_geographic_bounds = True

    # empty lists to hold data lines
    vx_list = []
    vy_list = []
    vz_list = []
    temp_list = []
    visc_list = []
    phase_list = []
    composition_list = []
    velocity_list = []

    # loop over caps
    for c in cap_list:
        c = '%02d' % int(c)

        # open the cap input file 
        inp = '%(datafile)s.cap%(c)s.%(time)s' % vars()
        inp_file = open(inp)
        if verbose: print (now(), 'xslice: open ', inp)

        # read the cap file into one big list of lines
        lines = inp_file.readlines()
        inp_file.close()
        if verbose: print (now(), 'xslice: close ', inp)

        # remove the header line
        lines = lines[1:]

        # composition file
        if 'comp' in field_list:
            print ("opening composition file")
            # form file name from model parameters
            opt = '%(datafile)s.opt%(c)s.%(time)s' % vars()
            # open the file 
            opt_file = open(opt)
            opt_lines = opt_file.readlines()
            # pop the header line
            opt_lines = opt_lines[1:]

        if verbose: print (now(), 'xslice: loop over file: ', inp)
        if verbose: print (now(), 'xslice: loop over j: 0,', nodey)
        if verbose: print (now(), 'xslice: loop over i: 0,', nodex)
        if verbose: print (now(), 'xslice: loop over k: 0,', nodez)
        # loop over the node numbers
        # NOTE loop order: Y , X , Z 
        for j in range(nodey):

            for i in range(nodex):

                # process only data at nodex == x 
                if not i == int(args['x']): 
                    continue # to next i 

                for k in range(nodez):

                    # compute line number
                    n = k + i*nodez + j*nodez*nodex

                    line = lines[n].rstrip()

                    # separate into cols: lat lon vx vy vz temp visc 
                    cols = line.split(' ')

                    # convert angles
                    fi = float(cols[1]) 
                    lon = fi * r2d

                    # non-dimensional radius
                    rad =      float(cols[2]) 
                    vx  =      float(cols[3])
                    vy  =      float(cols[4])
                    vz  =      float(cols[5])
                    temp =     float(cols[6])
                    visc =     math.log10(float(cols[7]))

                    if 'comp' in field_list:
                        opt_line = opt_lines[n].rstrip()


                    if check_geographic_bounds:
                        if r   < float(args['get_r_min']) \
                        or r   > float(args['get_r_max']) \
                        or lon < float(args['get_lon_min']) \
                        or lon > float(args['get_lon_max']) :
                            continue # to next k

                    # assemble data lists
                    temp_list.append( '%f %f %f' % (lon, rad, temp) )
                    visc_list.append( '%f %f %f' % (lon, rad, visc) )

                    if 'comp' in field_list:
                        composition_list.append( '%f %f %s' % (lon, rad,opt_line) )

                    # check for combined fields 
                    if field_list.count('velocity'):

                        # atan2 gives the proper quadrant 
                        # based on signs of vy and vz 
                        azimuth = math.degrees( math.atan2(vy,vz) )

                        # optionally, add the local lon to keep vectors 
                        # relative to local azimuth
                        if args.has_key('polar_angular_offset') :
                            # adjust azimuth for local position
                            azimuth += lon
                            # adjust azimuth for polar projection off set
                            a = float(args.get('polar_angular_offset'))
                            azimuth -= a

                        # apply scale 
                        length = (base_scale*math.hypot(vz,vy)) / velocity_scale
                        data = '%f %f %f %f' % (lon, rad, azimuth, length)
                        velocity_list.append( data ) 
                    # end of velocity calcs

                    # check for phase computations
                    if 'phase' in field_list:
                        phase = get_phase(phase_args, rad, temp) 
                        data = '%f %f %f' % (lon, rad, phase)
                        phase_list.append( data )

                # end of k loop over Z
            # end of i loop over X
        # end of j loop over Y
    # end of loop over caps


    # write out data files
    out_prefix = datafile + '.' \
                 + args['time'] + '.' \
                 + args['x'] + '.'

    if 'temp' in field_list:
        out_name = out_prefix + 'temp' + '.xyz'
        out_file = open(out_name, 'w')
        if verbose: print (now(), 'xslice: open,w ', out_name)
        for line in temp_list:
            out_file.write( '%s\n' % line )
        out_file.close()
        if verbose: print (now(), 'xslice: close: ', out_name)

    if 'visc' in field_list: 
        out_name = out_prefix + 'visc' + '.xyz'
        out_file = open(out_name, 'w')
        if verbose: print (now(), 'xslice: open,w ', out_name)
        for line in visc_list:
            out_file.write( '%s\n' % line )
        out_file.close()
        if verbose: print (now(), 'xslice: close: ', out_name)

    if 'phase' in field_list: 
        out_name = out_prefix + 'phase' + '.xyz'
        out_file = open(out_name, 'w')
        if verbose: print (now(), 'xslice: open,w ', out_name)
        for line in phase_list:
            out_file.write( '%s\n' % line )
        out_file.close()
        if verbose: print (now(), 'xslice: close: ', out_name)

    if 'comp' in field_list: 
        out_name = out_prefix + 'comp' + '.xyz'
        out_file = open(out_name, 'w')
        if verbose: print (now(), 'xslice: open,w ', out_name)
        for line in composition_list:
            out_file.write( '%s\n' % line )
        out_file.close()
        if verbose: print (now(), 'xslice: close: ', out_name)

    if 'velocity' in field_list: 
        out_name = out_prefix + 'velocity' + '.xyz'
        out_file = open(out_name, 'w')
        if verbose: print (now(), 'xslice: open,w ', out_name)
        p = 0
        q = 0 
        for line in velocity_list:
            q = q + 1
            # check for sub-sample
            if q % int(1/float(velocity_increment)) != 0:
                continue # to next line
            p = p + 1
            out_file.write( '%s\n' % line )
        out_file.close()
        if verbose: print (now(), 'xslice: total points read: ', q)
        if verbose: print (now(), 'xslice:    points written: ', p)
        if verbose: print (now(), 'xslice: close: ', out_name)


#=====================================================================
#=====================================================================
def yslice(args):
    '''Read a set of Citcoms cap files, write a set of .xyz files,
and return the set of file names as a comma separated list.
Output file name format: datafile.time.y.field_id
One output .xyz file will be produced for each field_id given.


Required arguments:
 args['citcoms_pid'] = name of pid file;
 args['time'] = model time steps [0,model end step];
 args['y'] = model node number in theta direction: [0, nodey - 1];
 args['field'] = comma separated list of field ids: 
 {temp, visc, velocity} (NOTE velocity parallel to slice plane);

Optional arguments: 
 args['velocity_scale'] = float value to scale velocity magnitude
NOTE: velocity data is produced in cm/yr multiplied by this factor.
 args['velocity_increment'] = float value, <= 1.0, to decimate data
 
 args['cap_list'] = used to explicitly set which caps to read and write

 args.get('get_r_min') = set limits on output region
 args.get('get_r_max') = set limits on output region
 args.get('get_lat_min') = set limits on output region
 args.get('get_lat_max') = set limits on output region

 args['polar_angular_offset'] = angle in degrees 
  Setting this value will cause the velocity vectors' azimuth value
  to be adjusted such that an azimuth of 0 will point 'locally up'.
  This is used to match the polar projection with GMT's PSXY.

Output arguments:
 none

Return value:
 list of output files written
'''

    # parse citcom model data
    citcoms = parse_citcoms_cfg_file( args['citcoms_pid'] ) 
    datafile = citcoms['datafile']
    nodex = int( citcoms['nodex'] ) 
    nodey = int( citcoms['nodey'] ) 
    nodez = int( citcoms['nodez'] ) 
    nproc_surf = int( citcoms['nproc_surf'] ) 

    # set cap list from nproc_surf
    cap_list = [0] 
    if nproc_surf == 12: 
        cap_list = range(12)

    # reset cap list if requested 
    if args.get('cap_list'):
        cap_list = args.get('cap_list')

    # parse time value
    time = args['time']

    # parse y value
    y = args['y']

    # parse field value, and settings
    value = args['field'] 
    # create a list with the value
    field_list = [ value ]
    # check for comma separated multiple values
    if value.count(','):
        # update list
        field_list = value.split(',')
    if verbose: print (now(), 'yslice: field_list ', field_list)

    # parse velocity scale
    velocity_scale = 1.0
    if args.get('velocity_scale'):
        velocity_scale = float(args.get('velocity_scale'))

    # convert velocities to cm/yr
    thermdiff = float( citcoms['thermdiff'] )
    layer_km = earth_radius
    base_scale = (thermdiff/(layer_km*1e3))*(100.0*3600.0*24.0*365.25)

    # parse velocity increment 
    velocity_increment = 1.0
    if args.get('velocity_increment'): 
        velocity_increment = float(args.get('velocity_increment'))

    # phase controls
    phase_args = args
    phase_args.update( citcoms )

    # parse geo bounds
    check_geographic_bounds = False
    if args.get('get_r_min') \
    or args.get('get_r_max') \
    or args.get('get_lon_min') \
    or args.get('get_lon_max') :
        check_geographic_bounds = True

    # empty lists to hold data lines
    vx_list = []
    vy_list = []
    vz_list = []
    temp_list = []
    visc_list = []
    phase_list = []
    composition_list =[] 
    velocity_list = []

    # loop over caps
    for c in cap_list:
        c = '%02d' % int(c)

        # open the cap input file 
        inp = '%(datafile)s.cap%(c)s.%(time)s' % vars()
        inp_file = open(inp)
        if verbose: print (now(), 'yslice: open ', inp)

        # read the cap file into one big list of lines
        lines = inp_file.readlines()
        inp_file.close()
        if verbose: print (now(), 'yslice: close ', inp)

        # remove the header line
        lines = lines[1:]

        # composition file
        if 'comp' in field_list:
            print ("opening composition file")
            # form file name from model parameters
            opt = '%(datafile)s.opt%(c)s.%(time)s' % vars()
            # open the file 
            opt_file = open(opt)
            opt_lines = opt_file.readlines()
            # pop the header line
            opt_lines = opt_lines[1:]

        if verbose: print (now(), 'yslice: loop over file: ', inp)
        if verbose: print (now(), 'yslice: loop over j: 0,', nodey)
        if verbose: print (now(), 'yslice: loop over i: 0,', nodex)
        if verbose: print (now(), 'yslice: loop over k: 0,', nodez)
        # loop over the node numbers
        # NOTE loop order: Y , X , Z 
        for j in range(nodey):

            # process only data at nodey == y 
            if not j == int(args['y']): 
                continue # to next j 

            for i in range(nodex):

                for k in range(nodez):

                    # compute line number
                    n = k + i*nodez + j*nodez*nodex

                    line = lines[n].rstrip()

                    # separate into cols: lat lon vx vy vz temp visc 
                    cols = line.split(' ')

                    # convert angles
                    theta = float(cols[0])
                    lat = 90 - (theta * r2d) 

                    # non-dimensional radius
                    rad =      float(cols[2]) 
                    vx  =      float(cols[3])
                    vy  =      float(cols[4])
                    vz  =      float(cols[5])
                    temp =     float(cols[6])
                    visc =     math.log10(float(cols[7]))

                    if 'comp' in field_list:
                        opt_line = opt_lines[n].rstrip()

                    if check_geographic_bounds:
                        if r   < float(args['get_r_min']) \
                        or r   > float(args['get_r_max']) \
                        or lat < float(args['get_lat_min']) \
                        or lat > float(args['get_lat_max']) :
                            continue # to next k

                    # assemble data lists
                    temp_list.append( '%f %f %f' % (lat, rad, temp) )
                    visc_list.append( '%f %f %f' % (lat, rad, visc) )

                    if 'comp' in field_list:
                        composition_list.append( '%f %f %s' % (lat,rad,opt_line) )

                    # check for combined fields 
                    if field_list.count('velocity'):

                        # atan2 gives the proper quadrant 
                        # based on signs of vy and vz 
                        azimuth = math.degrees( math.atan2(vx,vz) )

                        # optionally, add the local lon to keep vectors 
                        # relative to local azimuth
                        if args.has_key('polar_angular_offset') :
                            # adjust azimuth for local position
                            azimuth += lat
                            # adjust azimuth for polar projection off set
                            a = float(args.get('polar_angular_offset'))
                            azimuth -= a

                        # apply scale 
                        length = (base_scale*math.hypot(vz,vx)) / velocity_scale
                        data = '%f %f %f %f' % (lat, rad, azimuth, length)
                        velocity_list.append( data ) 
                    # end of velocity calcs

                    # check for phase computations
                    if 'phase' in field_list:
                        phase = get_phase(phase_args, rad, temp) 
                        data = '%f %f %f' % (lat, rad, phase)
                        phase_list.append( data )

                # end of k loop over Z
            # end of i loop over X
        # end of j loop over Y
    # end of loop over caps


    # write out data files
    out_prefix = datafile + '.' \
                 + args['time'] + '.' \
                 + args['y'] + '.'

    if 'temp' in field_list:
        out_name = out_prefix + 'temp' + '.xyz'
        out_file = open(out_name, 'w')
        if verbose: print (now(), 'yslice: open,w ', out_name)
        for line in temp_list:
            out_file.write( '%s\n' % line )
        out_file.close()
        if verbose: print (now(), 'yslice: close: ', out_name)

    if 'visc' in field_list: 
        out_name = out_prefix + 'visc' + '.xyz'
        out_file = open(out_name, 'w')
        if verbose: print (now(), 'yslice: open,w ', out_name)
        for line in visc_list:
            out_file.write( '%s\n' % line )
        out_file.close()
        if verbose: print (now(), 'yslice: close: ', out_name)

    if 'phase' in field_list: 
        out_name = out_prefix + 'phase' + '.xyz'
        out_file = open(out_name, 'w')
        if verbose: print (now(), 'yslice: open,w ', out_name)
        for line in phase_list:
            out_file.write( '%s\n' % line )
        out_file.close()
        if verbose: print (now(), 'yslice: close: ', out_name)

    if 'comp' in field_list: 
        out_name = out_prefix + 'comp' + '.xyz'
        out_file = open(out_name, 'w')
        if verbose: print (now(), 'yslice: open,w ', out_name)
        for line in composition_list:
            out_file.write( '%s\n' % line )
        out_file.close()
        if verbose: print (now(), 'yslice: close: ', out_name)

    if 'velocity' in field_list: 
        out_name = out_prefix + 'velocity' + '.xyz'
        out_file = open(out_name, 'w')
        if verbose: print (now(), 'yslice: open,w ', out_name)
        p = 0
        q = 0 
        for line in velocity_list:
            q = q + 1
            # check for sub-sample
            if q % int(1/float(velocity_increment)) != 0:
                continue # to next line
            p = p + 1
            out_file.write( '%s\n' % line )
        out_file.close()
        if verbose: print (now(), 'yslice: total points read: ', q)
        if verbose: print (now(), 'yslice:    points written: ', p)
        if verbose: print (now(), 'yslice: close: ', out_name)
#=====================================================================
#=====================================================================
def read_time_file(parameters_file):
    '''read a citcoms .time file; generate a list of tuples (step, age, runtime); where: step is model steps; age is in Ma; runtime in Myr.
Required arguments:
 parameters_file = name of pid file of model 

Optional arguments:
 none

Output arguments:
 none

Return value:
 list of time tuples 
''' 

    # establish a citcom parser
    parser = CitcomParser.Parser()
    parser.read(parameters_file)

    # get basic info about the case
    modelname = parser.getstr('datafile')
    start_age = parser.getfloat('start_age')
    thermdiff = parser.getfloat('thermdiff')
    output_format = parser.getstr('output_format')

    # compute temperature scale factor 
    # FIX where are these values from?
    layer_km = earth_radius
    scale_t = layer_km * 1e3 * layer_km * 1e3 / \
                  (thermdiff * 1.e6 * 365.25 * 24. * 3600.)
    
    # check for time file 
    # base name, local directory
    timefile = '%s.time' % modelname
    
    # check local dir first
    if not os.path.exists(timefile):

        # append sub dirs if local file not found:
        if output_format == 'hdf5':
            timefile = os.path.join('data', timefile)

        if output_format == 'ascii':
            timefile = os.path.join('data', '0', timefile)

        # re-check for file existance
        if not os.path.exists(timefile):
            raise (IOError, 'File not found: %(timefile)s') % vars()

    # open the file 
    try :
        time_file = open(timefile)
    except (IOError, (errno, strerror)):
        print ('I/O error(%s): %s') % (errno, strerror)
    except:
        print ('Unexpected error:', sys.exc_info()[0])
        raise

    # empty list
    tlist = []

    # read the file
    try: 
        lines = time_file.read().splitlines()
    finally:
        time_file.close()

    # loop over time file
    for line in lines:

        #step, t, dt, CPUtime, CPUdt = line.split() # FIX more general?
        step, t, dt, CPUtime, CPUdt = line.split(' ')
        step = int(step)

        # compute age in Ma
        age = start_age - ( float(t) * scale_t)
  
        # compute runtime in Myr
        runtime = float(t) * scale_t

        # fill the list
        tuple = (step, age, runtime)
        tlist.append(tuple)

    return tlist
#=====================================================================
#=====================================================================
def get_age_from_step(parameters_file, test_step):
    '''get age in Ma from model time step

Required arguments:
 parameters_file = pid file of model run
 test_step = model time step to convert to age

Optional arguments:
 none

Output arguments:
 none

Return value:
 age in Ma
''' 

    # get a list of tuples: (step, age, runtime) from the .time file
    list = read_time_file(parameters_file)

    # check bounds
    min = list[0][0]
    max = list[-1][0]
    if (int(test_step) > max):
        msg = 'step %s > than max step %s' % (test_step, max)
        raise (IndexError, msg)
    if (int(test_step) < min):
        msg = 'step %s < than min step %s' % (test_step, min)
        raise (IndexError, msg)

    # loop over tuples
    for (s,a,r) in list:
        if s == int(test_step):
            if verbose: 
                print (now(),'get_age_from_step: (s,a,r) =',(s,a,r))
            return a
    else:
        msg = 'step %s not found in .time file.' % test_step
        raise (IndexError, msg)    
#=====================================================================
#=====================================================================
def get_runtime_from_step(parameters_file, test_step):
    '''get runtime in Myr from model time step

Required arguments:
 parameters_file = pid file of model run
 test_step = model time step to convert to age

Optional arguments:
 none

Output arguments:
 none

Return value:
 run time in Myr 
'''

    # get a list of tuples: (step, age, runtime) from the .time file
    list = read_time_file(parameters_file)

    # check bounds
    min = list[0][0]
    max = list[-1][0]
    if (test_step > max):
        msg = 'step %s > than max step %s' % (test_step, max)
        raise (IndexError, msg)
    if (test_step < min):
        msg = 'step %s < than min step %s' % (test_step, min)
        raise (IndexError, msg)

    # loop over tuples
    for (s,a,r) in list:
        if s == int(test_step):
            if verbose:
                print (now(),'get_runtime_from_step: (s,a,r) =',(s,a,r))
            return r
    else:
        msg = 'step %(test_step)i not found in .time file.' % vars()
        raise (IndexError, msg)    

#=====================================================================
#=====================================================================
def get_step_from_age(parameters_file, test_age):
    '''locate and return the closest step for given age in Ma.

Required arguments:
 parameters_file = pid file of model run
 test_age = model age in Ma to convert to model time step

Optional arguments:
 none

Output arguments:
 none

Return value:
 integer model time step 

'''

    # get a list of tuples: (step, age, runtime) from the .time file
    list = read_time_file(parameters_file)

    # check bounds 
    max_age = list[0][1]
    min_age = list[-1][1]

    if (test_age > max_age):
        msg = 'age %f Ma > max age %f Ma' % (test_age, max_age)
        raise (IndexError, msg)

    if (test_age < min_age):
        msg = 'age %f Ma < min age %f Ma' % (test_age, min_age)
        raise (IndexError, msg)

    # shortcut for exact values
    for (s,a,r) in list:
        if a == test_age:
            if verbose: 
                print (now(),'get_step_from_age: (s,a,r) =',(s,a,r))
            return s

    # list comprehension to compute the difference between 
    # test_age and model ages
    delta_list = [a-test_age for (s,a,r) in list if a-test_age > 0 ]

    # get index of, and delta from value to test for previous value
    prev_i = len(delta_list) - 1
    prev_delta = delta_list[-1]

    # get index of, and delta from value to test for next value
    next_i = len(delta_list)
    next_delta = list[next_i][1] - test_age

    # index smallest delta and index into list of tuples
    if abs(prev_delta) < abs(next_delta) :
        i = prev_i
    else: 
        i = next_i

    # FIX 
    i = prev_i

    step = list[i][0]

    if verbose:
        print (now(), 'get_step_from_age:', 'prev=', list[prev_i], 'delta=', prev_delta, 'i=', prev_i)
        print (now(), 'get_step_from_age:', 'test=', test_age)
        print (now(), 'get_step_from_age:', 'next=', list[next_i], 'delta=', next_delta, 'i=', prev_i)

    return step
#=====================================================================
#=====================================================================
def get_step_from_runtime(parameters_file, test_runtime):
    '''locate and return the closest step for given runtime

Required arguments:
 parameters_file = pid file of model run
 test_runtime = model run time in Myr to convert to model time step

Optional arguments:
 none

Output arguments:
 none

Return value:
 integer model time step 
'''

    # get a list of tuples: (step, age, runtime) from the .time file
    list = read_time_file(parameters_file)

    # check bounds : NOTE indices
    max_runtime = list[-1][2]
    min_runtime = list[0][2]
    if (test_runtime > max_runtime):
        msg = 'runtime %f Myr > max runtime %f Myr' % (test_runtime, max_runtime)
        raise (IndexError, msg)
    if (test_runtime < min_runtime):
        msg = 'runtime %f Myr < min runtime %f Myr' % (test_runtime, min_runtime)
        raise (IndexError, msg)

    # shortcut for exact values
    for (s,a,r) in list:
        if float(r) == float(test_runtime):
            if verbose: 
                print (now(),'get_step_from_runtime: (s,a,r) =', (s,a,r))
            return s

    # list comprehension to compute the difference between 
    # test and values in list
    # NOTE : the if clause >
    delta_list = [r-test_runtime for (s,a,r) in list if r-test_runtime < 0 ]


    # get index of, and delta from value to test for previous value
    prev_i = len(delta_list) - 1
    prev_delta = delta_list[-1]

    # get index of, and delta from value to test for next value
    next_i = len(delta_list)
    next_delta = list[next_i][1] - test_runtime

    # index into main list of tuples
    i = 0
    if abs(prev_delta) < abs(next_delta) :
        i = prev_i
    else: 
        i = next_i

    step = list[i][0]
    if verbose:
        print (now(), 'get_step_from_runtime:', \
        'prev=', list[prev_i],'delta=', prev_delta, 'i=', prev_i)
        print (now(), 'get_step_from_runtime:', \
        'test=', test_runtime)
        print (now(), 'get_step_from_runtime:', \
        'next=', list[next_i],'delta=', next_delta,'i=', prev_i)
    return step
#=====================================================================
#=====================================================================
def get_time(pfile,time):
    '''get time in steps (0 to ..) from time in any string form: t in steps, tMa in age, tMyr in runtime

Required arguments:
 parameters_file = pid file of model run
 tim = string to convert to model time step

Optional arguments:
 none

Output arguments:
 none

Return value:
 model time steps
'''

    time = str( time )
    if time.endswith('Ma'):
        a = float(time[0:-2])
        t = get_step_from_age(pfile,a)
    elif time.endswith('Myr'):
        r = float(time[0:-3])
        t = get_step_from_runtime(pfile,r)
    else:
        t = int(time)
    return t
#=====================================================================
#=====================================================================
def read_citcoms_coord_0_file_into_zlist( dict ):
    '''read a global or regional case.coord.0 file; generate a list of tuples (level, z, depth); where: level ranges from 1 to nodez; z is the non-dimensional decimal fraction from model bottom (< 1.0) to surface (== 1.0); depth is in km.

Required arguments:
    dict['datafile']
    dict['datadir']
    dict['radius']

Optional arguments:
 none

Output arguments:
 none

Return value:
 list of vertical structure tuples 
'''
 
    # SAVE for later ... might need it ... 
    # identify regional or global case type (also called 'zone')
    # set header prefix string 
    #if parser.getint('nproc_surf') == 1:
    #    # regional
    #elif parser.getint('nproc_surf') == 12:
    #    # global
    #else:
    #    raise IOError, 'Unknown run type. \
    #                    For Regional set nproc_surf=1 or \
    #                    Global set nproc_surf=12 in \
    #                    %(parameters_file)' %vars()

    # get empty lists
    rlist = []
    zlist = []

    # variables from dictionary
    datafile = dict['datafile']
    datadir = dict['datadir']
    radius = float(dict['radius'])

    # correction if output by CitcomS processor specified in
    # pid file
    if datadir.endswith('%RANK'):
        # remove %RANK text
        datadir = datadir[:-len("%RANK")]
        # replace with 0 processor
        datadir = datadir + '0'

    # read the file into a list of lines
    file_name = '%(datadir)s/%(datafile)s.coord.0'  % vars()
    coord_file = open(file_name)
    try : 
        lines = coord_file.read().splitlines()
    finally:
        coord_file.close()

    # pop header line
    lines = lines[1:]

    z = 0 
    for line in lines:

        # parse the line
        (theta, phi, r) = line.split()
        r = float(r)

        if rlist.count(r) == 0:
            rlist.append(r)

            # compute depth in km
            depth = (1.0 - r) * radius*1E-3

            # fill the list
            tuple = (z, r, depth)
            zlist.append( tuple )
            z = z + 1
            if verbose: 
                print (now(), 'read_citcoms_coord_files: z, r, depth:', tuple)

        else: 
           continue # to next line
        
    return zlist

#=====================================================================
#=====================================================================
def read_citcoms_coor_file_into_zlist( dict ):
    '''read a global or regional coordinate file; generate a list of tuples (z, r, depth); where: z ranges from 0 to nodez-1; r is the non-dimensional radius, the decimal fraction from model bottom (< 1.0) to surface (== 1.0); depth is in km.

Required arguments:
 dict['coor']
 dict['coor_file']
 dict['nproc_surf']
 dict['radius']

Optional arguments:
 none

Output arguments:
 none

Return value:
 list of vertical structure tuples
'''
 
    coor = int(dict['coor'])
    nproc_surf = int(dict['nproc_surf'])
    coor_file = dict['coor_file']
    radius = float(dict['radius'])

    # if custom coor file not specified
    if coor == 0:
        # try reading the coord.0 file 
        return read_citcoms_coord_0_file_into_zlist( dict )

    # identify regional or global case type (also called 'zone')
    # set header prefix string 
    if nproc_surf == 1:
        # regional
        header_prefix = 'nsd= 3'
    elif nproc_surf == 12:
        # global
        header_prefix = 'nodez='
    else:
        raise (IOError, 'Unknown run type. \
                        For Regional set nproc_surf=1 or \
                        Global set nproc_surf=12 in \
                        %(parameters_file)') %vars()

    # read the file into a list of lines
    file = open(coor_file)
    try : 
        lines = file.read().splitlines()
    finally:
        file.close()

    # get an empty list of tuples (z, r, depth)
    zlist = []

    read = False

    for line in lines:

        # toggle read flag
        if line.lstrip().lstrip().startswith(header_prefix):
            read = not read 
            continue

        if read:

            # stop reading if another data block encountered 
            # delimited by lines starting with:
            #    ' nsd= N' where N is one of 1, 2, 3  
            # or 'nodeN' where N is one of 'x', 'y', 'z'
            if line.lstrip().lstrip().startswith('n'):
                read = not read 
                break;

            # parse the line
            (z, r) = line.split()
            z = int(z) - 1 
            r = float(r)
            
            # compute depth in km
            depth = (1.0 - r) * radius*1E-3

            # fill the list
            tuple = (z, r, depth)
            zlist.append( tuple )

            # WARN: uncomment for really verbose
            #if verbose: print now(), 'read_citcoms_coor_file: z, r, depth:', tuple

    return zlist
#=====================================================================
#=====================================================================
def get_depth_from_z(parameters_file, test_z):
    '''get depth in km from model z (0 .. nodez - 1)

Required arguments:
 parameters_file = pid file of model run 
 test_z = model node in z direction

Optional arguments:
 none

Output arguments:
 none

Return value:
 depth in km
'''

    test_z = int(test_z)

    # get a list of tuples: (z, r, depth) from coor file
    list = read_citcoms_coor_file_into_zlist(parameters_file)

    # check bounds
    min = list[0][0]
    max = list[-1][0]

    if (test_z > max):
        msg = 'z %s > max z %s' % (test_z, max)
        raise (IndexError, msg)
    if (test_z < min):
        msg = 'z %s < min z %s' % (test_z, min)
        raise (IndexError, msg)

    # loop over tuples
    for (z,r,d) in list:
        if test_z == int(z):
            return float(d)
    else:
        msg = 'z %s not found in coordinate file.' % z
        raise (IndexError, msg)
#=====================================================================
#=====================================================================
def get_r_from_z(parameters_file, test_z):
    '''get r, non-dimensional radius from model z (0 .. nodez - 1)

Required arguments:
 parameters_file = pid file of model run 
 test_z = model node in z direction

Optional arguments:
 none

Output arguments:
 none

Return value:
 non-dimentional radius 
'''

    test_z = int(test_z)

    # get a list of tuples: (z, r, depth) from coor file
    list = read_citcoms_coor_file_into_zlist(parameters_file)

    # check bounds
    min = list[0][0]
    max = list[-1][0]

    if (test_z > max):
        msg = 'z %s > max z %s' % (test_z, max)
        raise (IndexError, msg)
    if (test_z < min):
        msg = 'z %s < min z %s' % (test_z, min)
        raise (IndexError, msg)

    # loop over tuples
    for (z,r,d) in list:
        if test_z == int(z):
            return float(r)
    else:
        msg = 'z %s not found in coordinate file.' % z
        raise (IndexError, msg)
#=====================================================================
#=====================================================================
def get_z_from_depth(parameters_file, test_depth):
    '''locate and return closest z for given depth in km

Required arguments:
 parameters_file = pid file of model run 
 test_depth = depth in km

Optional arguments:
 none

Output arguments:
 none

Return value:
 model z node number
'''

    # get a list of tuples: (z, r, depth) from coor file
    zlist = read_citcoms_coor_file_into_zlist(parameters_file)

    # check bounds
    max = zlist[0][2]
    min = zlist[-1][2]
    if (test_depth > max):
        msg = 'depth %s km > max depth %s km' % (test_depth, max)
        raise (IndexError, msg)
    if (test_depth < min):
        msg = 'depth %s km < min depth %s km' % (test_depth, min)
        raise (IndexError, msg)

    # shortcut for exact values
    for (z,r,d) in zlist:
        if d == float(test_depth): 
            if verbose: 
                print (now(),'get_z_from_depth: (z,r,d) =', (z,r,d))
            return z
        
    # list comprehension to compute the difference between
    # test_depth and model depths for all deltas in list
    # NOTE : the if clause > 
    delta_list = [ d - test_depth for (z,r,d) in zlist if d - test_depth > 0 ] 

    # get index of, and delta from depth value previous to test_depth
    prev_i = len(delta_list) - 1
    prev_delta = delta_list[-1]

    # get index of, and delta from depth value next to test_depth
    next_i = len(delta_list)
    next_delta = zlist[next_i][1] - test_depth

    # find smallest delta
    if abs(prev_delta) < abs(next_delta) :
        i = prev_i
    else:
        i = next_i

    z = zlist[i][0]
    if verbose: 
        print (now(),'get_z_from_depth:',\
        'prev=', zlist[prev_i],'delta=',prev_delta,'i=',prev_i)
        print (now(),'get_z_from_depth:',\
        'test=', test_depth)
        print (now(),'get_z_from_depth:',\
        'next=', zlist[next_i],'delta=',next_delta,'i=',next_i)
        print (now(),'get_z_from_depth:',\
        'index=', i, 'z=', zlist[i][0])
    return z
#=====================================================================
#=====================================================================
def get_z_from_r(parameters_file, test_r):
    '''locate and return closest z for given r in non-dimensional radius coordinates (0.55 to 1.0)

Required arguments:
 parameters_file = pid file of model run 
 test_r = non-dimensional radius value

Optional arguments:
 none

Output arguments:
 none

Return value:
 model z node number
'''

    zlist = read_citcoms_coor_file_into_zlist(parameters_file)

    # check bounds
    max = zlist[-1][1]
    min = zlist[0][1]
    if (test_r > max):
        msg = 'r value %s > than max r %s' % (test_r, max)
        raise (IndexError, msg)
    if (test_r < min):
        msg = 'r value %s < than min r %s' % (test_r, min)
        raise (IndexError, msg)

    # shortcut for exact values
    for (z,r,d) in zlist:
        if float(z) == float(test_r): 
            if verbose: 
                print (now(),'get_z_from_r: (z,r,d) = ', (z,r,d))
            return z
        
    # list comprehension to compute the difference between
    # test_r and model z for all deltas in list
    # NOTE : the if clause < 
    delta_list = [ z - test_r for (l,z,d) in zlist if z - test_r < 0 ] 

    # get index of, and delta from z value previous to test_r
    prev_i = len(delta_list) - 1
    prev_delta = delta_list[-1]

    # get index of, and delta from z value next to test_r
    next_i = len(delta_list)
    next_delta = zlist[next_i][1] - test_r

    # find smallest delta
    if abs(prev_delta) < abs(next_delta) :
        i = prev_i
    else:
        i = next_i

    z = zlist[i][0]
    if verbose: 
        print (now(),'get_z_from_r:',\
        'prev=', zlist[prev_i],'delta=',prev_delta,'i=',prev_i)
        print (now(),'get_z_from_r:',\
        'test=', test_r)
        print (now(),'get_z_from_r:',\
        'next=', zlist[next_i],'delta=',next_delta,'i=',next_i)
        print (now(),'get_z_from_r:',\
        'index=', i, 'level=', zlist[i][0])
    return z
#=====================================================================
#=====================================================================
def get_phase(args, radius, t):
    '''From a dictionary containing CitcomS model parameters, a non-dimensional radius value and a temperature value, compute total phase.

Required arguments:
 args.get('radius_outer') = outer radius of the model, in non-dimensional units
 radius = radius value to find phase for 
 t = model time step 

Optional arguments:
 none

Output arguments:
 none

Return value:
 total phase value 
'''

    # outer radius of the model, in non-dimensional units
    radius_outer = args.get('radius_outer')
    if radius_outer == None: 
        msg = "radius_outer must be set in pid file"
        raise (ValueError, msg)

    # convert strings to numbers
    radius_outer = float ( radius_outer )

    #if verbose: print now(), "get_phase: radius_outer =", radius_outer
    #if verbose: print now(), "get_phase: radius =", radius
    #if verbose: print now(), "get_phase: t =", t

    # initial total phase value
    total_phase = 0.0

    # loop over phase ids
    for id in ['410', '670', 'cmb']:

        # get the phase parameters for this id
        width = args.get('width' + id)
        transT = args.get('transT' + id)
        clapeyron = args.get('clapeyron' + id)
        Ra_phase = args.get('Ra_' + id)

        # skip missing sections
        if width == None \
        or transT == None \
        or Ra_phase == None \
        or clapeyron == None :
            continue # to next id

        # convert strings to numbers
        width = float ( width )
        transT = float ( transT )
        clapeyron = float ( clapeyron )
        Ra_phase = float ( Ra_phase )

        # skip ids with Ra_phase == zero
        if int(Ra_phase) == 0:
            #if verbose: print now(), "get_phase: skip this id = ", id
            continue # to next id

        # NOTE name change.
        if id == '670':  
            z_level = args.get('z_' + 'lmantle')
        else:
            z_level = args.get('z_' + id)

        # convert strings to numbers
        z_level = float ( z_level )

        # compute pressure differential and phase

        pressure = (radius_outer - radius) - z_level - \
                       clapeyron * (t - transT)

        phase = 0.5 * ( 1.0 + math.tanh( pressure / width ) )

        total_phase += phase

        #if verbose: print now(), "get_phase: Ra_phase = ", Ra_phase
        #if verbose: print now(), "get_phase: z_level =", z_level
        #if verbose: print now(), "get_phase: clapeyron =", clapeyron
        #if verbose: print now(), "get_phase: transT =", transT
        #if verbose: print now(), "get_phase: width =", width
        #if verbose: print now(), 'get_phase: pressure =', pressure
        #if verbose: print now(), 'get_phase: phase =', phase
        # end of loop over ids

    #if verbose: print now(), 'get_phase: total_phase =', total_phase

    return total_phase

#=====================================================================
#=====================================================================
#=====================================================================
#=====================================================================
#=====================================================================
def test_time_functions():
    '''test of time functions'''

    # parse cmd line
    file = sys.argv[1]
    test_step = int( sys.argv[2] )
    test_age = float( sys.argv[3] )
    test_runtime = float( sys.argv[4] )
    
    age = get_age_from_step( file, test_step )
    print (now(), 'test_time_functions: test_step = %(test_step)i --> age = %(age)f') % vars()

    runtime = get_runtime_from_step( file, test_step )
    print (now(), 'test_time_functions: test_step = %(test_step)i --> runtime = %(runtime)f') % vars()

    step = get_step_from_age( file, test_age )
    print (now(), 'test_time_functions: test_age = %(test_age)f --> step = %(step)i') % vars()

    step = get_step_from_runtime( file, test_runtime )
    print (now(), 'test_time_functions: test_runtime = %(test_runtime)f --> step = %(step)i') % vars()

#=====================================================================
#=====================================================================
def test_depth_functions():
    '''test of depth functions'''

    pid_file = sys.argv[1]

    dict = parse_citcoms_cfg_file(pid_file)

    # test depth from z for all z, 0 to nodez - 1 
    for z in range( int(dict['nodez']) ):
        d = get_depth_from_z(pid_file, z)
        print (now(), 'test_depth_functions: get_depth_from_z(pid, z) z = %(z)s; d = %(d)s') % vars()

    print (' ') 
    print (' ') 

    # test depth from z for all z, 0 to nodez - 1 
    for z in range( int(dict['nodez']) ):
        r = get_r_from_z(pid_file, z)
        print (now(), 'test_depth_functions: get_r_from_z(pid, z) z = %(z)s; r = %(r)s') % vars()

    print (' ') 
    print (' ') 

    # test z from depth for some values, 0 to 2000
    d_list = range(0, 2000, 200)
    for d in d_list:
        z = get_z_from_depth(pid_file, d) 
        print (now(), 'test_depth_functions: get_z_from_depth(pid_file, d) z = %(z)s; d = %(d)s') % vars()
        print (' ') 

    print (' ') 
    print (' ') 

    # test z from r for some values from 0.6 to 1.0 
    r_list = [i * 0.1 for i in range(6, 11)]
    for r in r_list:
        z = get_z_from_r(pid_file, r) 
        print (now(), 'test_depth_functions: get_z_from_r(pid_file, d) z = %(z)s; r = %(r)s') % vars()
        print (' ') 


#=====================================================================
#=====================================================================
def test_cap_polygons( argv ):
    pid_file = sys.argv[1]
    time = sys.argv[2]
    d = {} 
    d['citcoms_pid'] = pid_file
    d['time'] = time 
    d['cap_list'] = [0]
    d['cap_list'] = range(12)

    write_cap_polygons( d )
#=====================================================================
#=====================================================================
def test_slice( argv ):
    '''test slice functions'''
    global verbose
    verbose = True 

    pid_file = sys.argv[1]
    time = sys.argv[2]
    n = sys.argv[3]
    field = sys.argv[4]

    dict = {} 

    # set data read parameters:
    dict['citcoms_pid'] = pid_file
    dict['time'] = time 
    dict['z'] = n
    dict['y'] = n
    dict['x'] = n
    dict['field'] = field

    dict['velocity_scale'] = 1.0
    dict['velocity_increment'] = 1.0

    # sub-test for limited cap list 
    # dict['cap_list'] = [4] 

    #dict['slice_type'] = 'zslice'

    # control slice type
    dict['slice_type'] = 'xslice'

    slice( dict )

#=====================================================================
#=====================================================================

def Sfile( dict ):
    '''Makes a topography xy file 
Required arguments:
    dict containing the following values
    1. dict['refvisc']
    2. dict['thermdiff']
    3. dict['density']
    4. dict['gravacc']
    5. dict['datafile']
    6. dict['nodex']
    6. dict['nodey']

Must have a suface file in the working directory
which is output from CitcomS.
This file has the string *surface* in it. There should
only be one in your working directory.

Adds the following file to the dictionary:
    dict['surf_xy_file'] = surffile

Optional arguments:
 none

Output arguments:
 none

Return value:
 none
    '''
    # Make scale factor for conversion from non-dimensional vertical stress
    # You need to dimensionalise then convert to a factor to find h.
    refvisc=float(dict['refvisc'])
    therm_diff = float(dict['thermdiff'])
    radius=float(6371.0)
    density=float(dict['density'])
    densitya=float(0.0)
    densityw=float(1000.0)
    delrho=density-densityw
    gravacc=float(dict['gravacc'])
    stress_scale=1/((1/(delrho*gravacc)) * ((refvisc*therm_diff)/(radius*radius)))
    # Define scale for scaling dynamic topography
    scale_s=stress_scale
    # Define output surface xyz file
    surffile=dict['datafile']+'_surf.xys' 
    r2d = 180.0/math.pi
    # Find the surface file in your current working directory
    # sfile=commands.getoutput('ls *surface*') # pyhton 2 version
    sfile=subprocess.getoutput('ls *surface*') # python 3 version, replaced commands with subprocess O.F.B 2023
    # make arrays of node numbers
    nodex=float(dict['nodex'])
    nodey=float(dict['nodey'])
    nodes_per_proc=nodex*nodey
    inputd=open(sfile,"r")
    output=open(surffile,"w")
    line=inputd.readline()
    # Write an xyz file for dynamic topography with scaled data.
    for n in range(nodes_per_proc):
        line=inputd.readline()
        longs=(float(line.split(' ')[1]))*r2d
        lat=90 - ((float(line.split(' ')[0])) * r2d)
        surf=float(line.split(' ')[2])
        surf=scale_s*surf				
        output.write('%f %f %f \n' % (longs,lat,surf))
    inputd.close()
    output.close()
    dict['surf_xy_file'] = surffile
    print ('made surface file...')
    return


#=====================================================================

#=====================================================================
#=====================================================================
def map_global_proc_to_proc_in_z_dir( args ):
    '''Create a list that maps global proc to proc in z dir 

Required arguments:
 args['nproc_surf'] = no. of caps
 args['nprocx']     = no. of processors in x direction
 args['nprocy']     = no. of processors in y direction
 args['nprocz']     = no. of processors in z direction

Optional arguments:
 none

Return:
 none    
'''
    # get variables from dictionary
    nproc_surf = int(args['nproc_surf'])
    nprocx = int(args['nprocx'])
    nprocy = int(args['nprocy'])
    nprocz = int(args['nprocz'])

    # make list
    proc_z_level = []

    for n in range( nproc_surf ):
        for j in range( nprocy ):
            for i in range( nprocx ):
                for k in range( nprocz ):
                    proc_z_level.append(k)

    # update dictionary
    args['map_list'] = proc_z_level

    return

#=====================================================================
#=====================================================================
def map_track_file_data_to_proc_node_number_list( args ):
    '''Map the track file data to a global node number list

Required arguments:
 args['nx']          = no. of processors in the x direction
 args['ny']          = no. of processors in the y direction
 args['nz']          = no. of processors in the z direction
 args['total_procs'] = total no. of processors
 args['track_files'] = list of track files arranged by z node number

Optional arguments:
 none

Return value:
 none
'''

    # get variables from dictionary
    nx = int(args['nx'])
    ny = int(args['ny'])
    nz = int(args['nz'])
    track_files = args['track_files']
    total_procs = int(args['total_procs'])

    # initialise the list
    list = []

    # set the length of the list (so that cells can be
    # directly referenced by their list index)
    list = [0 for i in range( nx*ny*nz )]

    for k in range( nz ):
        input = open(track_files[k]) # open track file at depth
        lines = input.readlines()
        input.close
        linenum = 0 # counter for line number

        for j in range( ny ):
            for i in range( nx ):
                node = k + i*nz + j*nz*nx # N.B. k and NOT (k+1)
                temp = float(lines[linenum].split()[2])
                list[node] = temp
                linenum += 1 # increment counter

    # update dictionary
    args['temp_list'] = list

    return
#=====================================================================
#=====================================================================
def test( argv ):
    '''self test function'''

    global verbose
    verbose = True
    # test time 
    #test_time_functions()

    # test depth functions
    test_depth_functions()

    # test slicing 
    # test_slice( argv )

    # test cap polygon 
    # test_cap_polygons( argv )
#=====================================================================
if __name__ == "__main__":
    import Core_Citcom

    if len( sys.argv ) > 1:
        # process sys.arv as file name(s) for testing 
        test( sys.argv )
    else:
        # print module documentaiton and exit
        help(Core_Citcom)
#=====================================================================
#=====================================================================

#!/usr/bin/env python
#=====================================================================
#                      Python Scripts for CitComS 
#         Preprocessing, Data Assimilation, and Postprocessing
#                  ---------------------------------
#
#                    AUTHORS: Mark Turner, Dan Bower
#
#             (c) California Institute of Technology 2007
#                        ALL RIGHTS RESERVED
#=====================================================================
'''A set of general purpose functions for working with Files.

This module has functions to read and parse various file types associated 
with CitComS, GPlates, the plot_page system, and other data sets.

Most functions take a dictionary of arguments to pass and adjust paramters.
'''
#=====================================================================
#=====================================================================
# imports
import getopt, os, string, sys, math, time
import datetime
import pprint

import CitcomParser 
import Core_Citcom
import Core_GMT
from Core_Util import now
#=====================================================================
#=====================================================================
# global variables
verbose = False
#=====================================================================
#=====================================================================

#=====================================================================
def read_slab_db_into_dictionary( filename ):
    '''read a slab_db file into a dictionary and return it'''

    if verbose:
        print now(), 'read_slab_db_into_dictionary:', filename

    file = open(filename)
    try : 
        # read into a list of lines
        lines = file.read().splitlines()
    finally:
        file.close()

    # main slab dictionary to construct ; a dictionary of dictionaries
    slab_dictionary = {}

    # the key for the main slab_dictionary
    slab_name = ''    

    # the value for the main slab_dictionary
    sub_dict = {} 

    # read
    for line in lines:

        # skip blank lines
        if string.strip(line) == '':
            continue # to next line 

        if not line.startswith('#'):

            # save this line as a slab name, a key to the main slab dict
            slab_name = line
            continue # to next line 

        if line.startswith('#'):

            # get the sub-dictionary data by splitting on # 
            element_array = line.split( '#' )

            # loop over the data elements in the line 
            for e in element_array:

                # skip empty elements from extra #'s at start & end of line
                if e == '' :
                    continue; # to next element

                if e.find('=') == -1:
                    # 
                    # NOTE: first element of slab_db are not need; 
                    # do not add them to the dictionary , but save this code for reference
                    #

                    # # this element is the plate polygon name
                    # key = 'plate_polygon'
                    # val = e
                    # # remove white space
                    # val = val.strip()
                    # val = val.rstrip()
                    # # populate the sub-dictionary
                    # sub_dict[key] = val
                    continue; # to next element

                else :
                    [key,val] = e.split('=')
                    # remove white space
                    key = key.strip()
                    key = key.rstrip()

                    # NOTE: these elements of slab_db are not need; 
                    if key == 'poly_start_age': continue; # to next element 
                    if key == 'poly_end_age': continue; # to next element 

                    val = val.strip()
                    val = val.rstrip()

                    # populate the sub-dictionary
                    sub_dict[key] = val

            # deep copy this sub_dict to the main slab dictionary
            slab_dictionary[slab_name] = dict(sub_dict) 

            # clear out old values of this working sub_dict
            sub_dict = {}

            continue # to next line 

        else:
           print "WARNING: line not understood as key or values?:"
           print line


    if verbose:
      print now(), 'read_slab_db_into_dictionary: slab_dictionary='

      pprint.PrettyPrinter(indent=2).pprint(slab_dictionary)

    return slab_dictionary
#=====================================================================
#=====================================================================
def parse_control_file( filename ):
    '''Parse a key=value style control file, and return a dictionary of those settings; Each section of the control file will be a separate sub-dictionary.'''

    # empty settings dict
    settings = {}

    # section flag data
    section = None
    num_sections = 0

    file = open(filename)

    try : 
        # read into a list of lines
        lines = file.read().splitlines()
    finally:
        file.close()

    if verbose:
        print now(), 'parse_control_file: read: %(filename)s' % vars()

    # read
    for line in lines:

        # skip blank lines
        if string.strip(line) == '':
            continue # to next line in control file

        # skip lines starting with '#'
        if line.startswith('#'):
            continue # to next line in control file

        # re-set the section flag on end of section lines
        if line.startswith('[END') : 
            section = None
            if verbose:
                print now(), 'parse_control_file: END section = ', section

            continue # to next line in control file

        # process section boundary into new key 
        if line.startswith('[') :
            line = line.rstrip()
            line = line.replace('[','')
            line = line.replace(']','')
            section = line
            num_sections += 1
            if verbose:
                print now(), 'parse_control_file: section = ', section

            # establish an empty map for this section
            settings[section] = {}

            continue # to next line in control file

        # key/value pairs can be seperated by whitespaces
        opt = line.rstrip()
        opt = opt.replace(' ','') # remove any spaces padding

        # isolate the key and value
        (key,val) = string.split(opt, '=')

        if verbose:
            print now(), 'parse_control_file: key , val = ',key,',',val

        if (section) :
            # this key val pair is part of a figure spec
            settings[section][key] = val

        else : 
            # this key val pair is a general setting
            settings[key] = val

    if verbose:
        print now(), 'parse_control_file: settings ='
        pprint.PrettyPrinter(indent=2).pprint(settings)
        print now(), 'parse_control_file: settings ='

    return settings
#=====================================================================
#=====================================================================
def gplates_velocity_to_xyz(args):
    '''merge CitComS coordinte data with GPlates generated velocity data, and return the filename'''

    # constant of nature 
    r2d = 180.0/math.pi

    # determin how many velocity files to process
    cap_num = 0 
    if args['zone'] == 'global':
        cap_num = 12
    else:
        cap_num = 1
    
    iage = args['iage']
    if not args.get('overlay_gplates_velocity_increment'):
        msg = 'gplates_velocity_to_xyz: ERROR: Parameter "overlay_gplates_velocity_increment" is reqired to read gplates velocity files' % vars()
        raise ValueError, msg

    if not args.get('overlay_gplates_velocity_vector_scale') :
        msg = 'gplates_velocity_to_xyz: ERROR: Parameter "overlay_gplates_velocity_vector_scale" is reqired to read gplates velocity files' % vars()
        raise ValueError, msg

    scale = args['overlay_gplates_velocity_vector_scale']
    increment = args['overlay_gplates_velocity_increment']

    # final output file
    out = 'gplates_velocity.%s.xyz' % (iage)
    out_file = open(out, 'w')
    if verbose: print now(), 'cap_to_xyVelocity: open,w ', out

    # set initial path to find gplates data 
    path = args.get('overlay_gplates_velocity_path')
    if path == None:
        path = './'

    model = args.get('overlay_gplates_velocity_coord')

    # loop over caps
    for c in range (cap_num):
        cor = '%(model)s.coord.%(c)s' % vars()
        tmp1 = '%(model)s.coord.%(c)s.tmp1' % vars()
        tmp2 = '%(model)s.coord.%(c)s.tmp2' % vars()

        # build gplates file name
        prefix = args.get('overlay_gplates_velocity_prefix')
        if prefix.startswith('/'):
            path = ''
        f = prefix + '%(iage)s.%(c)s' % vars() 
        gpv = os.path.join(path, f )

        if verbose: print now(), 'gplates_velocity_to_xyz: cor =', cor
        if verbose: print now(), 'gplates_velocity_to_xyz: tmp1 =', tmp1
        if verbose: print now(), 'gplates_velocity_to_xyz: tmp2 =', tmp2
        if verbose: print now(), 'gplates_velocity_to_xyz: gpv =', gpv

        # check for files
        if not os.path.exists(gpv):
            msg = 'gplates_velocity_to_xyz: ERROR: file missing: %(gpv)s' % vars()
            raise ValueError, msg
        
        # check for files
        if not os.path.exists(cor):
            msg = 'gplates_velocity_to_xyz: ERROR: file missing: %(cor)s' % vars()
            raise ValueError, msg

        # build temporay file 
        cmd = 'grep "1.000000e+00" %(cor)s | cut -d" " -f1,2,3 > %(tmp1)s' % vars()
        if verbose: print now(), 'gplates_velocity_to_xyz: cmd =', cmd
        os.system(cmd)

        cmd = 'paste -d" " %(tmp1)s %(gpv)s > %(tmp2)s' % vars()
        if verbose: print now(), 'gplates_velocity_to_xyz: cmd =', cmd
        os.system(cmd)

        # open the input file 
        # read the inp file into one big list of lines
        inp_file = open(tmp2)
        if verbose: print now(), 'gplates_velocity_to_xyz: read ', tmp2
        lines = inp_file.readlines()
        inp_file.close()
        if verbose: print now(), 'gplates_velocity_to_xyz: close ', tmp2
        if verbose: print now(), 'gplates_velocity_to_xyz: loop start '
        p = 0
        q = 0
        for line in lines:
            q = q + 1
            # sub-sample data 
            if q % int(1/float(increment)) == 0: 
                   p = p + 1
                   #if verbose: print now(), 'cap_to_xyVelocity: q========', q
                   line = line.rstrip()
                   # separate into cols: lat lon rad vx vy 
                   cols = line.split()

                   # convert angles
                   lat = 90 - float(cols[0]) * r2d
                   lon =      float(cols[1]) * r2d
                   rad =      float(cols[2]) 
                   vx  =      float(cols[3])
                   vy  =      float(cols[4])
      
                   if vy <= 0.0: 
                       if vx < 0.0:
                           a=math.degrees(math.atan(vy/vx))
                           azimuth=360.0-a
                       elif vx > 0.0:
                           a=math.degrees(math.atan(-vy/vx))
                           azimuth=180.0+a
                       else:
                           azimuth=270.0
                   elif vy > 0.0:
                       if vx < 0.0:
                           a=math.degrees(math.atan(-vy/vx))
                           azimuth=a
                       elif vx > 0.0:
                           a=math.degrees(math.atan(vy/vx))
                           azimuth=180.0-a
                       else:
                           azimuth=90.0

                   # scale to get inches on page
                   length = math.hypot(vx,vy) / float(scale)

                   # write data
                   out_file.write( '%f %f %f %f\n' % (lon, lat, azimuth, length) )
            # end of if q mod incre == 0 increment sub-sampling
        # end of loop over lines
        if verbose: print now(), 'gplates_velocity_to_xyz: total =', q
        if verbose: print now(), 'gplates_velocity_to_xyz: used =', p
        if verbose: print now(), 'gplates_velocity_to_xyz: loop end '

        cmd = 'rm -v %(tmp1)s %(tmp2)s' % vars()
        if verbose: print now(), 'gplates_velocity_to_xyz: cmd =', cmd
        os.system(cmd)

    # end of loop over caps

    # close file 
    out_file.close()
    if verbose: print now(), 'gplates_velocity_to_xyz: close: ', out

    return out

#=====================================================================
#=====================================================================
#=====================================================================
def export_citcom_xy_nodal_coord_by_proc( args ):
    '''Export CitcomS xy (lon,lat) nodal coord data by processor.

Required arguments:
 args['datafile']    = CitcomS file prefix
 args['nx']          = no. of x nodes per processor
 args['ny']          = no. of y nodes per processor
 args['nz']          = no. of z nodes per processor
 args['total_procs'] = total no. of processors
 args['proc']        = processor no. 

Optional arguments:
 none

Return value:
 none
'''

# Read in CitcomS coor file proc by proc and then
# create new temporary xy file that contains the nodal coordinates
# of just one level that can be used by GMT.

    # Set handy variable
    r2d = 180.0/math.pi

    # Pull out required variables from dictionary

    nx = int(args['nx']) 
    ny = int(args['ny'])
    nz = int(args['nz'])

    datafile=args['datafile']

    total_procs = int(args['total_procs'])
    proc = int(args['proc']) # processor number

    # required for looping procedure
    num = nx * ny

    # maybe we should make the location of the
    # coordinate files a user defined variable in
    # the input cfg file.
    # e.g. coord_files = /this/is/an/example
    ifile = 'Coord/%s.coord.%d' % (datafile,proc) # input coord file
    input = open(ifile)

    # read file and close original
    lines = input.read().split('\n')
    input.close()
    lines = lines[1:] # pop header

    # output file name for nodel coordinate (NC) file
    ofile = 'NC.%d.xy' % (proc)
    output = open(ofile,'w')

    # extract and write out data
    for q in range ( num ):
        theta = float(lines[nz*q].split(' ')[0])
        phi = float(lines[nz*q].split(' ')[1])
        lat = 90. - float(theta) * r2d
        lon = float(phi) * r2d
        if lon < 0:
            lon += 360

        output.write( '%g %g\n' % (lon,lat) )

    output.close()

    # update dictionary
    args['xy'] = ofile

    return

#=====================================================================
#=====================================================================
#=====================================================================
def grdtrack_grid_files_by_depth( args ):
    '''Track through a list of grid files in the x-y plane for each
 citcom z node (i.e. depth) and write out track files.

Required arguments:
 args['nz']          = no. of z nodes per processor
 args['total_procs'] = total no. of processors
 args['proc']        = processor no.
 args['grid_files']  = list of grid files (ordered by citcom z node)
 args['map_list']    = list of mapping from global proc to proc in z dir
                       see Core_Citcom.map_global_proc_to_proc_in_z_dir

Optional arguments:
 none

Return value:
 none
'''
    # get variables from dictionary
    nz = int(args['nz'])
    total_procs = int(args['total_procs'])
    proc = int(args['proc'])
    map_list = args['map_list']

    track_files = [] # list to store names of track files

    for k in range ( nz ):
        # kk is the z node number
        kk = map_list[proc]*(nz-1)+k
        # the grd file given by argument 'G' will be a function
        # of z node i.e. kk.  Currently this variable is hard-coded
        # for testing ONLY!
        args['G'] = args['grid_files'][kk] # select grd file
        track_file = Core_GMT.grdtrack( args ) # track file name
        track_files.append( track_file )

    # update dictionary
    args['track_files'] = track_files

    return

#=====================================================================
#=====================================================================
def write_citcom_velo ( args ):
    '''parse temperature distribution in the correct ordering
to the citcom 'velo' files.

Required arguments:
 args['nx']          = no. of x nodes per processor
 args['ny']          = no. of y nodes per processor
 args['nz']          = no. of z nodes per processor
 args['total_procs'] = total no. of processors
 args['proc']        = processor no.
 args['temp_list']   = list of temperatures ordered by proc
                       node number

Optional arguments:
 none

Return value:
 none
'''

   # get variables from dictionary
    nx = int(args['nx'])
    ny = int(args['ny'])
    nz = int(args['nz'])
    total_procs = int(args['total_procs'])
    proc = int(args['proc'])
    total_nodes_proc = nx*ny*nz
    temp_list = args['temp_list']

    # velo file
    velofile = 'IC/test.velo.%(proc)s.0' % vars()
    vfile = open(velofile,'w')

    izero = 0
    ione = 1

    # velo file header
    vfile.write('%(izero)d %(total_nodes_proc)d %(izero)e\n' % vars() )
    vfile.write('  %(ione)d %(total_nodes_proc)d\n' % vars() )

    for j in range( ny ):
        for i in range( nx ):
            for k in range( nz ):
                node = k + i*nz + j*nz*nx
                temp = temp_list[node]
                vfile.write('%(izero)e %(izero)e %(izero)e %(temp)e\n' % vars () )


    vfile.close()

    return

#=====================================================================
#=====================================================================
def export_to_citcom_velo( args ):
    '''Export grid files to citcom velo files (by proc).

Required arguments:
 args['nodex']       = no. of nodes in x dir
 args['nodey']       = no. of nodes in y dir
 args['nodez']       = no. of nodes in z dir
 args['nprocx']      = no. of processors in x dir
 args['nprocy']      = no. of processors in y dir
 args['nprocz']      = no. of processors in z dir
 args['total_procs'] = total no. of processors
 args['nx']          = no. of x nodes per processor
 args['ny']          = no. of y nodes per processor
 args['nz']          = no. of z nodes per processor
 args['grid_files']  = a list of grid files that correspond by index to
                       the citcom z node.  For example, grid_files[0]
                       is the grid file associated with the 'first'
                       citcom z node in the radial direction

Optional arguments:
 none

Return value:
 none
'''

    # get variables from dictionary
    total_procs = int(args['total_procs'])

    # this function also sets args['map_list'] with the mapping
    # list data
    if verbose: print now(), 'Core_File: export_to_citcom_velo: call Core_Citcom.map_global_proc_to_proc_in_z_dir'
    Core_Citcom.map_global_proc_to_proc_in_z_dir( args )

    # master loop over all processors to export data
    for proc in range( total_procs ):

        # processor number required by subsequent functions
        args['proc'] = proc

        # this function also sets args['xy'] = output file
        # args['xy'] is required by subsequent functions
        if verbose: print now(), 'Core_File: export_to_citcom_velo: call Core_File.export_citcom_xy_nodal_coord_by_proc'
        export_citcom_xy_nodal_coord_by_proc( args )

        # this function also sets args['track_files'] = output files
        # args['track_files'] is required by subsequent functions
        if verbose: print now(), 'Core_File: export_to_citcom_velo: call Core_File.grdtrack_grid_files_by_depth'
        grdtrack_grid_files_by_depth( args )

        # this function also sets args['temp_list'] = output temperature
        # args['temp_list'] is required by export_to_citcom_velo
        if verbose: print now(), 'Core_File: export_to_citcom_velo: call Core_Citcom.map_track_file_data_to_proc_node_number_list'
        Core_Citcom.map_track_file_data_to_proc_node_number_list( args )

        # export velo files 
        if verbose: print now(), 'Core_File: export_to_citcom_velo: call Core_File.write_citcom_velo'
        write_citcom_velo( args )

        # export mat files?

    return

#=====================================================================
#=====================================================================
def test( argv ):
    '''self test'''
    global verbose
    verbose = True 

    # args = {} 
    # print now(), 'test: parse_citcoms_cfg_file(%(pid_file)s)'
    # args = parse_citcoms_cfg_file(pid_file)

# Sample test code
    # dict = {} 
    # print now(), 'test: parse_control_file(%(file)s)'
    # dict = parse_control_file(file)
#=====================================================================
if __name__ == "__main__":
    import Core_File

    if len( sys.argv ) > 1:
        # process sys.arv as file names for testing 
        test( sys.argv )
    else:
        # print module documentaiton and exit
        help(Core_File)
#=====================================================================
#=====================================================================

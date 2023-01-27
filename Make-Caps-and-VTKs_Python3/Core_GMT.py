#!/usr/bin/python
#=====================================================================
#                      Python Scripts for CitComS 
#         Preprocessing, Data Assimilation, and Postprocessing
#                  ---------------------------------
#             (c) California Institute of Technology 2008
#                        ALL RIGHTS RESERVED
#=====================================================================
''' A set of functions to build and execute GMT commands.

This module has functions to call GMT commands with versitle options.

Most functions take a dictionary of arguments to pass and adjust paramters.

'''
#=====================================================================
#=====================================================================
# imports
import os, string, sys, math, time, datetime, pprint, random, commands

import Core_Util
from Core_Util import now

#=====================================================================
# Global variables
verbose = False
#=====================================================================
#=====================================================================
# funtions 
#=====================================================================
#=====================================================================
def mask_grid_by_nan_grid(args):
    '''Mask out an input grid by a second grid with NAN values

Required arguments:
 args['grid'] = input grid to mask
 args[out'] = output grid from mask
 args['mask_grid'] = grid file containing areas of NANs and data
    
Optional arguments:

Output arguments:
 none

Return value:
 none
'''

    # set up file names
    grid = args['grid'] # the input grid to process
    out = args['out'] # the output grid to generate
    mask = args['mask_grid'] # the masking grid with NANs

    file = mask[ mask.rfind('/')+1 : ]
    isnan = 'tmp_isnan_copy_of' + file
    clip = 'tmp_clip_from_' + isnan

    # Use grdmath to generate a grid with only 1 or 0 values
    cmd = 'grdmath -V %(mask)s ISNAN = %(isnan)s' % vars()
    print now(), 'mask_grid_by_nan_grid: cmd =\n', cmd
    os.system(cmd)

    # Use grdclip to change the 0's to NAN's
    cmd = 'grdclip -V %(isnan)s -Sb1/NAN -G%(clip)s' % vars()
    print now(), 'mask_grid_by_nan_grid: cmd =\n', cmd
    os.system(cmd)

    # Use grdmath to 'mask' the input grid by the NAN grid
    cmd = 'grdmath -V %(grid)s %(clip)s OR = %(out)s' % vars()
    print now(), 'mask_grid_by_nan_grid: cmd =\n', cmd
    os.system(cmd)

    # clean up
    cmd = 'rm -rf %(isnan)s %(clip)s' % vars()
    print now(), 'mask_grid_by_nan_grid: cmd =\n', cmd
    os.system(cmd)

#=====================================================================
#=====================================================================
def mask_grid_by_xy_file(args):
    '''Mask out grid data from polygonal areas.

Required arguments:
 args['grid'] = input grid to mask
 args[out'] = output grid from mask
 args['xy'] = xy file containing closed polygon outline data to use as mask

Optional arguments:

Output arguments:
 none

Return value:
 none
'''
    # xy file with the paths to mask from 
    xy = args['xy'] 

    # the source grid to mask ; it will set -R 
    grid = args['grid'] 

    # the final grid with values inside the mask pasths
    out = args['out']
    
    # other local vars for the GMT calls 

    # the tmp grid from the mask 
    tmp_grid_from_mask = 'tmp_grid_from_mask.grd' 

    # 1.get basic info about it into args
    #grdinfo( args )

    # 2. generate the grid with the masking
    # NOTE: -R from source grid 
    # NOTE: -N out/edge/in
    cmd = 'grdmask %(xy)s -G%(tmp_grid_from_mask)s -R%(grid)s -V -m -NNaN/1/1 ' % vars()
    # will lat,lon vs. lon,lat ever go away?
    if args.has_key(':') :
      toggle = args[':'] 
      cmd += '-%(toggle)s' % vars()
    print now(), 'mask_grid_by_xy_file: cmd =\n', cmd
    os.system(cmd)

    # Apply the xy mask file to src grid using grdmath. 
    cmd = 'grdmath -V %(grid)s %(tmp_grid_from_mask)s MUL = %(out)s' % vars()
    print now(), 'mask_grid_by_xy_file: cmd =\n', cmd
    os.system(cmd)

    # clean up
    cmd = 'rm -rf %(tmp_grid_from_mask)s' % vars()
    print now(), 'mask_grid_by_xy_file: cmd =\n', cmd
    os.system(cmd)

#=====================================================================
def start_page(args):
    '''Establish a .ps or .eps file by calling gmtset and pstext.

Required arguments:
 args['ps'] = name of file to create 
    
Optional arguments:
 args['eps']
 args['page_title'] = text string to print at top of page
 args['page_orientation'] = 'landscape' or 'portrait' (default)

Output arguments
 args['K'] = '-K'
 args['O'] = ''
 args['R'] = '0/11.0/0/8.5' or '0/8.5/0/11.0 -P ' (default)

Return value:
 none
'''
    # remove any existing page 
    ps = args['ps']
    cmd = 'rm -rf %(ps)s\n' % vars()
    if verbose: print now(), 'start_page: cmd=\n', cmd
    os.system(cmd)

    # remove any existing .gmtdefaults 
    ps = args['ps']
    cmd = 'rm -rf .gmtdefaults4\n' % vars()
    if verbose: print now(), 'start_page: cmd=\n', cmd
    os.system(cmd)

    # check for eps
    cmd = ''
    if args.get('eps'):
        cmd = 'gmtset PAPER_MEDIA letter+\n'
    else:
        cmd = 'gmtset PAPER_MEDIA letter\n'

    cmd += 'gmtset MEASURE_UNIT inch\n'

    for key in args.keys():
        if key.startswith('gmtset '):
            cmd += '%s %s\n' % ( key, args[key] )

    cmd += '\n'

    if verbose: print now(), 'start_page: cmd=\n', cmd
    os.system(cmd)
    
    text = '  '
    if args.get('page_title'):
        text = args['page_title'] 

    # set up page text 
    if args.get('page_orientation') == 'landscape':
        # landscape
        args['R'] = '0/11.0/0/8.5'
        args['text'] = '0 7.0 18 0 0 BL %(text)s' % vars()
    else:
        # portrait
        args['R'] = '0/8.5/0/11.0 -P '
        args['text'] = '0 9.0 18 0 0 BL %(text)s' % vars()

    # keep postscript open option
    args['K'] = ' '
    pstext(args)

#=====================================================================
#=====================================================================
def end_page(args):
    '''End a ps file by calling pstext 

Required arguments:
 args['ps'] = name of file to create 
    
Optional arguments:
 args['page_number'] = text string to print at bottom of page;

Output arguments
 args['R'] = '0/11.0/0/8.5' or '0/8.5/0/11.0 -P ' (default)

Return value:
 none
'''

    if verbose: print now(), 'end_page:'

    text = ' ' 
    if args.get('page_number'):
        text = args['page_number'] 
    else:
        text = '_' 

    # set up page footer 
    if args.get('page_orientation') == 'landscape':
        # landscape
        args['text'] = '0 0 10 0 0 BL %(text)s' % vars()
        args['R'] = '0/11.0/0/8.5'

    else:
        # portrait
        args['text'] = '0 0 10 0 0 BL %(text)s' % vars()
        args['R'] = '0/8.5/0/11.0 -P '

    # close postscript 
    args['O'] = ' ' # ensure overlay 
    if args.has_key('K'):
        del args['K'] # continuation no longer required
    pstext(args)

#=====================================================================
#=====================================================================
def psbasemap( args ):
    '''call psbasemap

Required arguments:
 args['J']
 args['B']
 args['ps']
 args['R']

Optional arguments:
 args['O']
 args['K']
 args['X']
 args['Y']
 args['V']

Output arguments:

Return value:
 none
'''

    B = args['B']
    J = args['J']
    R = args['R']
    ps = args['ps']

    cmd = 'psbasemap -B%(B)s -J%(J)s -R%(R)s ' % vars()

    if verbose: cmd += '-V '
    if args.get('X'): cmd += '-X%s ' % args.get('X') 
    if args.get('Y'): cmd += '-Y%s ' % args.get('Y') 
    if args.get('K'): cmd += '-K%s ' % args.get('K') 
    if args.get('O'): cmd += '-O%s ' % args.get('O') 

    cmd += '>> %(ps)s' % vars()

    if verbose: print now(), 'psbasemap: cmd=\n', cmd
    os.system(cmd)

#=====================================================================
#=====================================================================
def pstext(args):
    '''call pstext

Required arguments:
 args['text'] = text to print
 args['R'] = page region
 args['ps'] = file to write text into 
    
Optional arguments:
 args['O'] = continuation flag
 args['K'] = continuation flag
 args['G'] = pen weight code
 args['W'] = Paint a rectangle beneath the text string. Set color [Default is no fill]

Output arguments:
 none 

Return value:
 none 
'''

    text = args['text']
    R = args['R']
    ps = args['ps']

    cmd = 'pstext -N -Jx1.0 -R%(R)s ' % vars()

    if verbose: cmd += '-V '
    if args.get('G'): cmd += '-G%s ' % args.get('G') 
    if args.get('W'): cmd += '-W%s ' % args.get('W') 
    if args.get('K'): cmd += '-K%s ' % args.get('K') 
    if args.get('O'): cmd += '-O%s ' % args.get('O') 

    cmd += '<<END >> %(ps)s\n' % vars()

    cmd +='%(text)s\n' % vars()
    cmd +='END'
    if verbose: print now(), 'pstext: cmd =\n', cmd
    os.system(cmd)

#=====================================================================
#=====================================================================
def grdtrack( args ):
    '''Call grdtrack

Required arguments:
 args['xy'] = xy file
 args['G'] = grd file

Optional arguments:
 none

Return value:
 track_file = 'G.track'
'''

    xy = args['xy']
    G = args['G']

    track_file = "%(G)s.track" % vars()

    cmd = 'grdtrack %(xy)s -G%(G)s > %(track_file)s' % vars()

    if verbose: cmd += ' -V'

    if verbose: print now(), 'grdtrack: cmd= ', cmd
    os.system(cmd)

    return track_file

#=====================================================================
#=====================================================================
def get_region( args ):
    '''Call minmax on an xyz file to set the region, -R.

Required arguments:
 args['xyz_file'] = name of xyz file to parse;

Output arguments:
 args['R'] = w/e/s/n of xyz file;

Output arguments:
 none 

Return value:
 none 
'''

    xyz = args['xyz_file']

    # query rotated xyz for -R values
    cmd = 'minmax -I0.01/0.01 %(xyz)s' % vars()
    if verbose: print now(), 'get_region: cmd =\n', cmd
    pipe = os.popen(cmd)
    line = pipe.readline()
    R = line.strip()
    pipe.close()
    # strip off preceeding '-R' 
    R = R[2:]
    if verbose: print now(), 'get_region: R=', R
    args['R'] = R
#=====================================================================
#=====================================================================
def xyz2grd( args ): 
    '''backward compatiblity wrapper to call make_grid

NOTE: at somepoint this wrapper will be filled with 
actual call to GMT's xyz2grd
'''

    f = make_grid( args )
    return f

#=====================================================================
def make_grid( args ) :
    '''Create a grid file using GMT calls.

Required arguments:
 args['xyz_file'] = xyz file to make grid from;
 args['grid_increment']  sets -I value
 args['grid_min'] sets -Ll value
 args['grid_max'] sets -Lu value
 args['grid_tension'] sets -T value

Optional arguments:
 args['headers'] sets -H value
 args['file_in_latlon'] sets -: flag

Output arguments:
 None

Return value:
 grid file = 'xyz_file.grd'
'''

    # file parameters
    xyz_file = args['xyz_file']
    mean_file="%s_median.xyz" % xyz_file
    grid_file = "%s.grd" % xyz_file

    # geography
    R = args['R']

    # grid parameters
    grid_increment = args['grid_increment']
    grid_min = args['grid_min']
    grid_max = args['grid_max']
    grid_tension = args['grid_tension']

    # smooth the data into a mean file 
    cmd = 'blockmedian %(xyz_file)s ' % vars()
    if verbose: cmd += '-V ' 
    if args.get('file_headers'): 
        H = args.get('file_headers')
        cmd += '-H%(H)s ' % vars()
    if args.get('file_in_latlon'): cmd += '-: ' % vars()
    cmd += '-I%(grid_increment)s -R%(R)s > %(mean_file)s' % vars()
    if verbose: print now(), 'make_grid: cmd =\n', cmd
    os.system(cmd)

    # check for file 
    if not os.path.exists( mean_file ):
        msg = 'blockmean may have failed: file not found: %(mean_file)s' % vars()
        raise ValueError, msg

    # create the grid
    cmd = 'surface %(mean_file)s ' % vars()
    if verbose: cmd += '-V '
    if args.get('file_headers'): 
        H = args.get('file_headers')
        cmd += '-H%(H)s ' % vars()
    if args.get('file_in_latlon'): cmd += '-: ' % vars()
    if grid_min != 'none' or grid_max != 'none':
        cmd += '-Ll%(grid_min)s -Lu%(grid_max)s ' % vars()
    cmd += '-G%(grid_file)s -I%(grid_increment)s -R%(R)s -T%(grid_tension)s -V ' % vars()
    if verbose: print now(), 'make_grid: cmd =\n', cmd
    os.system(cmd)

    # check for file 
    if not os.path.exists( grid_file ):
        msg = 'surface may have failed: file not found: %(mean_file)s' % vars()
        raise ValueError, msg

    # clean up
    cmd = "rm -rf %(mean_file)s" % vars()
    if verbose: print now(), 'make_grid: cmd =\n', cmd
    os.system(cmd)

    return grid_file
#=====================================================================
#=====================================================================
def grdinfo(args):
    '''Get basic info about a grid and update the args dictionary.

Required arguments:
 args['grid'] = name of grid file to parse 

Optional arguments: 
 None

Output argugments:
 args['west']
 args['east']
 args['south']
 args['north']
 args['z0']
 args['z1']
 args['dx']
 args['dy']
 args['nx']
 args['ny']
 args['x0']
 args['y0']
 args['x1']
 args['y1']
 args['n_nan']
 
 args['R'] = w/e/s/n
'''        

# NOTE:
#        -C produces info on a single line using the format:
#           file w e s n z0 z1 dx dy nx ny [x0 y0 x1 y1] [n_nan]
#                1 2 3 4 5  6  7  8  9  10 11  12 13 14  15
#
    grid = args['grid']

    cmd = 'grdinfo '
    if verbose: cmd += '-V '
    cmd += '-C -M %(grid)s' % vars() 
    if verbose: print now(), 'grdinfo: cmd =\n', cmd
    pipe = os.popen(cmd)
    line = pipe.readline()
    pipe.close()

    list = line.split()
    args['west'] = list[1]
    args['east'] = list[2]
    args['south'] = list[3]
    args['north'] = list[4]
    args['z0'] = list[5]
    args['z1'] = list[6]
    args['dx'] = list[7]
    args['dy'] = list[8]
    args['nx'] = list[9]
    args['ny'] = list[10]
    args['x0'] = list[11]
    args['y0'] = list[12]
    args['x1'] = list[13]
    args['y1'] = list[14]
    args['n_nan'] = list[15]

    if verbose:
        print 'west =', list[1]
        print 'east =', list[2]
        print 'south =', list[3]
        print 'north =', list[4]
        print 'z0 =', list[5]
        print 'z1 =', list[6]
        print 'dx =', list[7]
        print 'dy =', list[8]
        print 'nx =', list[9]
        print 'ny =', list[10]
        print 'x0 =', list[11]
        print 'y0 =', list[12]
        print 'x1 =', list[13]
        print 'y1 =', list[14]
        print 'n_nan =', list[15]
        print 'R  =', str( '/'.join(list[1:5]) )

    # combine raw grid data into plot parameters
    args['R'] = '/'.join(list[1:5]) 
#=====================================================================
#=====================================================================
def makecpt(args):
    '''make a cpt file, either from scratch or a given cpt file.

Required arguments:
 none
    
Optional arguments:
 args['cpt_file'] = if supplied, use value as the basis for the new cpt 

 args['C'] = if supplied, use value as the basis for the new cpt, 
  else use 'polar'

 args['cpt_z0'] or args['z0'] or args['grid_min'] = set min value 
 args['cpt_z1'] or args['z1'] or args['grid_max'] = set max value 
 args['cpt_dz'] or args['dz'] or args['grid_increment'] = set delta

 args['I'] = Invert the scale 
 args['Z'] = Make a continuous scale

Output arguments
 args['C'] = new cpt file 

Return value:
 none

'''
    # first, check for existing file name; 
    if args.has_key('cpt_file'):
        # update args
        args['C'] = cpt_file
        return
    else:
        # generate a random name
        cpt_file = 'tmp.' + str(os.getpid()) + '.cpt'

    # check for explicit base pallet
    C = args.get('C')
    if not C:
        C = 'polar'

    # for each setting, z0, z1, dz, 
    # first check for explicit value
    # next, check for value from grid settings
    # finally, check for value from grid info 
    use_T = True

    if args.has_key('cpt_z0'):
        z0 = args.get('cpt_z0')
    elif args.has_key('grid_min'):
        z0 = args.get('grid_min')
    elif args.has_key('z0'):
        z0 = args['z0']
    else:
        use_T = False
        
    if args.has_key('cpt_z1'): 
        z1 = args.get('cpt_z1')
    elif args.has_key('grid_max'):
        z1 = args.get('grid_max')
    elif args.has_key('z1'):
        z1 = args['z1']
    else:
        use_T = False

    if args.has_key('cpt_dz'):
        dz = args.get('cpt_dz')
    elif args.has_key('grid_increment'):
        dz = args.get('grid_increment')
    elif args.has_key('dz'):
        dz = args['dz']
    else:
        use_T = False

    # cmd and required arguments
    cmd = 'makecpt -M '
    if use_T: 
        cmd += '-T%(z0)s/%(z1)s/%(dz)s ' % vars()
    cmd += '-C%(C)s ' % vars()

    # optional arguments
    if args.get('V') or verbose: 
        cmd += '-V ' % vars()
    if args.get('I') or args.get('cpt_invert'):
        cmd += '-I ' % vars()
    if args.get('M'):
        cmd += '-M ' % vars()
    if args.get('Z') or args.get('cpt_cont'): 
        cmd += '-Z ' % vars()

    # cmd output
    cmd += '> %(cpt_file)s' % vars()
    if verbose: print now(), 'makecpt: cmd =\n', cmd
    os.system(cmd)

    # update args
    args['C'] = cpt_file
#=====================================================================
#=====================================================================
def psscale( args ):
    '''plot a scale bar

Required arguments:
 args['ps'] = psfile to plot scale bar on
 args['C'] = cpt file
    
Optional arguments:
 args['B'] = explicit scale lable, and boundary info 
 args['scale_label'] = scale label , required if B not set 
 args['D'] 

Output arguments:
 none

Return value:
 none
'''

    # required args
    ps = args.get('ps')

    # Check for explict values
    D = args.get('D')
    if not D:
        D = '0/0/2/0.25h'

    B = args.get('B')
    if not B:
        B = 'a1.0f0.5:%s:' % args['scale_label']

    # cmd and required args
    C = args.get('C')

    cmd = 'psscale '
    cmd += '-D%(D)s -B%(B)s -C%(C)s ' % vars()

    # optional args
    if args.get('K'): cmd += '-K%s ' % args.get('K') 
    if args.get('O'): cmd += '-O%s ' % args.get('O') 
    if args.get('Q'): cmd += '-Q%s ' % args.get('Q')

    if args.get('V') or verbose: 
        cmd += '-V ' % vars()

    # cmd output 
    cmd += ' >> %(ps)s' % vars()
    if verbose: print now(), 'psscale: cmd =\n', cmd
    os.system(cmd)
#=====================================================================
#=====================================================================
def grdimage(args):
    '''Create grayshaded or colored image from a 2-D netCDF grid file

Required arguments:
 args['grid'] = grid file to plot
 args['ps'] = ps file to draw into
 args['J'] = projection code
 args['C'] = cpt file 
    
Optional arguments:
 args['R'] = Region
 args['B'] = Boundary info
 args['X'] = X position on page
 args['Y'] = Y position on page
 args['K'] = GMT Keep ps file open
 args['O'] = GMT Overlay mode
 args['V'] = GMT verbose code
 args['xy'] = plot an xy overlay

Output arguments:
 none

Return value:
 none
'''


    grid = args['grid']
    ps = args['ps']
   
    # command with required arguments
    cmd = 'grdimage %(grid)s ' % vars()

    cmd += '-J%s ' % args['J']
    cmd += '-C%s ' % args['C']

    if verbose: cmd += '-V ' 
    if args.get('V'): cmd += '-V%s ' % args.get('V') 
    if args.get('B'): cmd += '-B%s ' % args.get('B') 
    if args.get('R'): cmd += '-R%s ' % args.get('R') 
    if args.get('X'): cmd += '-X%s ' % args.get('X') 
    if args.get('Y'): cmd += '-Y%s ' % args.get('Y') 
    if args.get('K'): cmd += '-K%s ' % args.get('K') 
    if args.get('O'): cmd += '-O%s ' % args.get('O') 

    cmd += ' >> %(ps)s' % vars()
    if verbose: print now(), 'grdimage: cmd =\n', cmd
    os.system(cmd)
#=====================================================================
#=====================================================================
def grdcontour(args):
    '''Contouring of 2-D gridded data sets

Required arguments:
 args['grid'] = grid file to plot
 args['ps'] = ps file to draw into
 args['J'] = projection code
 args['C'] = cpt file 
    
Optional arguments:
 args['R'] = Region
 args['B'] = Boundary info
 args['X'] = X position on page
 args['Y'] = Y position on page
 args['K'] = GMT Keep ps file open
 args['O'] = GMT Overlay mode
 args['V'] = GMT verbose code
 args['A'] = Annotation interval code
 args['W'] = Pen weight code

Output arguments:
 none

Return value:
 none
'''

    grid = args['grid']
    ps = args['ps']
   
    # command with required arguments
    cmd = 'grdcontour %(grid)s ' % vars()

    cmd += '-J%s ' % args['J']
    cmd += '-C%s ' % args['C']

    if verbose: cmd += '-V ' 
    if args.get('V'): cmd += '-V%s ' % args.get('V') 
    if args.get('B'): cmd += '-B%s ' % args.get('B') 
    if args.get('R'): cmd += '-R%s ' % args.get('R') 
    if args.get('X'): cmd += '-X%s ' % args.get('X') 
    if args.get('Y'): cmd += '-Y%s ' % args.get('Y') 
    if args.get('K'): cmd += '-K%s ' % args.get('K') 
    if args.get('O'): cmd += '-O%s ' % args.get('O') 
    if args.get('W'): cmd += '-W%s ' % args.get('W')
    if args.get('A'): cmd += '-A%s ' % args.get('A')

    cmd += ' >> %(ps)s' % vars()
    if verbose: print now(), 'grdcontour: cmd =\n', cmd
    os.system(cmd)
#=====================================================================
#=====================================================================
def psxy(args):
    '''Plot lines, polygons, and symbols on maps

Requried arguments:
 args['ps'] = ps file for output
 args['xy'] = xy file to plot
 args['R'] = region
 args['J'] = projection

Optional arguments:
 args['G'] = Color code 
 args['S'] = symbol codes
 args['M'] = use a single blank space for multi segment file 
 args[':'] = use a single blank space for data in lat lon

Output arguments:
 none

Return value:
 none
'''
    xy = args['xy']
    ps = args['ps']
   
    # command with required arguments
    cmd = 'psxy %(xy)s ' % vars()
    cmd += '-J%s ' % args['J']
    cmd += '-R%s ' % args['R']

    # optional args
    if verbose: cmd += '-V ' 
    if args.get('B'): cmd += '-B%s ' % args.get('B') 
    if args.get('G'): cmd += '-G%s ' % args.get('G') 
    if args.get('K'): cmd += '-K%s ' % args.get('K') 
    if args.get('M'): cmd += '-m%s ' % args.get('M') 
    if args.get('m'): cmd += '-m%s ' % args.get('m') 
    if args.get('O'): cmd += '-O%s ' % args.get('O') 
    if args.get('S'): cmd += '-S%s ' % args.get('S') 
    if args.get('W'): cmd += '-W%s ' % args.get('W') 
    if args.get('X'): cmd += '-X%s ' % args.get('X') 
    if args.get('Y'): cmd += '-Y%s ' % args.get('Y')
    if args.get('N'): cmd += '-N '  
    if args.get(':'): cmd += '-:%s ' % args.get(':') 

    cmd += ' >> %(ps)s' % vars()
    if verbose: print now(), 'psxy: cmd =\n', cmd
    os.system(cmd)
#=====================================================================
#=====================================================================
def grdrotater( args ):
    '''Rotate a grid using a finite rotation

Required arguments:
 args['grid'] = input grid to rotate
 args['out'] = output grid from rotation
 args['lon'] = Euler Pole lon
 args['lat'] = Euler Pole lat
 args['angle'] = Euler Pole rotation angle
   
Optional arguments:
 none

Output arguments:
 none

Return value:
 none
'''

    in_grid = args['grid']
    out_grid = args['out']

    lon = args['lon']
    lat = args['lat']
    angle = args['angle']
   
    # outline of projected region
    #outline  = out_grid + '_projected_outline_in.xy'

    # command with required arguments
    cmd = 'grdrotater %(in_grid)s ' % vars()

    cmd += '-T%(lon)s/%(lat)s/%(angle)s ' % vars()
    cmd += '-G%(out_grid)s ' % vars()
    cmd += '-N ' % vars()
    cmd += '-Q ' % vars()

    #if args.has_key('R'):
    #  R = args['R']
    #  cmd += '-R%(R)s ' % vars()

    if verbose: 
        cmd += '-V ' % vars()
        print now(), 'grdrotater: cmd =\n', cmd
    # SAVEME: outline 
    # cmd += '> %(outline)s' % vars()
    os.system(cmd)
#=====================================================================
#=====================================================================
def pslegend(args):
    '''create postscript for legend

Required arguments:
 args['ps'] = psfile to draw on 
 args['pslegend_file'] = text file containing pslegend code lines
 args['pslegend_D'] or args['D'] = position info for legend
    
Optional arguments:
 args['F'] = add a frame around the legend
 args['page_orientation'] = set page type, default is portrait

Output arguments:
 none

Return value:
 none
'''

    ps = args['ps']

    txt = args['pslegend_file']

    if args.get('pslegend_D'):
        D = args.get('pslegend_D') 
    elif args.get('D'): 
        D = args.get('D') 
    else:
        msg = "must set args['pslegend_file'] and args['pslegend_D'] or args['D'] before calling pslegend" % vars()
        raise ValueError, msg

    cmd = 'pslegend ' % vars()

    if verbose: cmd += '-V '

    if args.get('K'): cmd += '-K%s ' % args.get('K') 
    if args.get('O'): cmd += '-O%s ' % args.get('O') 

    # optionally draw a frame around the legend 
    if args.get('F'): cmd += '-F%s ' % args.get('F') 

    # set region from page orientation
    R = '0/8.5/0/11.0'

    if args.get('page_orientation') == 'portrait':
        R = '0/8.5/0/11.0'
    elif args.get('page_orientation') == 'landscape':
        R = '0/11/0/8.5'

    cmd += '-R%(R)s -Jx1i -Gwhite -D%(D)s "%(txt)s" >> %(ps)s' % vars()
    if verbose: print now(), 'pslegend: cmd =\n', cmd
    os.system(cmd)
#=====================================================================
#=====================================================================
def write_pslegend_file(args):
    '''Write a .txt file with figure information in pslegend format.

Required arguments:
 args['ps'] = used to name the pslegend_file 
 args['figure_id'] = used to name the pslegend_file 

Optional arguments:
 args['overlay_legend_info']
 args['overlay_velocity_vector_legend']
 args['overlay_velocity_vector_scale']
 args['overlay_gplates_velocity_vector_legend']
 args['overlay_gplates_velocity_vector_scale'] 
 args['overlay_platepolygons'] 

 etc. ... see python code for more options...

Output arguments:
 args['pslegend_file'] = file_name of generated text 
 args['pslegend_D'] = automatically calculated height, depending on number and kind of optional section given above

Return value:
 none
'''

    # ps and figure_id are used to associate legend file with ps file
    ps = args['ps'] 
    figure_id = args['figure_id'] 

    datafile = args['datafile']

    t = args.get('time')
    l = args.get('level')

    J = args['J']
    R = args['R']

    # compute the height of the legend as items are added 
    h = 0.0

    # assemble the full pslegend text string 
    txt = '' % vars()

    # Run info 
    if args.get('overlay_legend_info'): 
        txt += '''# info line
L 9 Helvetica L %(figure_id)s: %(datafile)s, time=%(t)s, level=%(l)s
G 0.05i
D 0 1p,black,
''' % vars()
        h += 0.20

    #
    # Velocity vector
    #
    if args.get('overlay_velocity_vector_legend'): 
        s = args.get('overlay_velocity_vector_scale') or 1.0
        txt += '''# velocity vector
S 0.5 v 1.0/0.015/0.06/0.05 0/0/0 1,black, 1.3i Velocity scale: %(s)s cm/yr 
G 0.05i
D 0 1p,black,
''' % vars()
        h += 0.25

    # GPlates Velocity vector
    if args.get('overlay_gplates_velocity_vector_legend'): 
        s = args.get('overlay_gplates_velocity_vector_scale') or 1.0
        txt += '''# velocity vector
S 0.5 v 1.0/0.015/0.06/0.05 0/0/0 1,black, 1.3i Velocity scale: %(s)s cm per year 
G 0.05i
D 0 1p,black,
''' % vars()
        h += 0.25


    # Feature data 
    if args.get('overlay_feature_legend'):
        # start legend section
        txt += '''# feature types
''' % vars()

        # all line data
        if args.get('overlay_lines'):
            txt += '''# line data
S 0.1i - 0.15i - 3p,black, 0.3i Feature Data
D 0 1p,black,
''' % vars()
            h += 0.20

        # all plate polygon data
        if args.get('overlay_platepolygons'):
            txt += '''# feature types
S 0.1i - 0.15i - 3p,green, 0.3i Plate Polygons 
D 0 1p,black,
''' % vars()
            h += 0.20

        # boundary types
        if \
        args.get('overlay_subduction_boundaries_sL') or \
        args.get('overlay_subduction_boundaries_sR') or \
        args.get('overlay_ridge_transform_boundaries'): 
            txt += '''# plate boundary types
N 3
V 0 1p,black,
S 0.1i - 0.15i - 3p,black, 0.3i sR
S 0.1i - 0.15i - 3p,black, 0.3i sL
S 0.1i - 0.15i - 3p,green, 0.3i Ridge/Trans.
V 0 1p,black,
D 0 1p,black,
N 1
''' % vars()
            h += 0.20

        # individual plate polygon outlines
        if args.get('overlay_platepolygon_boundaries') :
            # parse string for plate ids
            list = args['overlay_platepolygon_boundaries'].split(',')
            txt += '''# plate polygons
N 3
V 0 1p,black,
''' % vars()
            i = 0
            for prefix in list:
                i += 1

                # locate the color
                pid = prefix[ prefix.rfind('.')+1: ]
                color = PlatesColorTable.get_color(pid)

                # COLOR Correct 
                if color.upper() == 'YELLOW':
                    color = 'dark' + color 

                # add a line for each plate polygon
                txt += '''S 0.1i - 0.15i - 3p,%(color)s, 0.3i Plate Id %(pid)s
''' % vars()
                # add a new height delta every 3 plate id's 
                if i % 4 == 0: 
                   h += 0.20

            # finish out plate polygon section
            txt += '''V 0 1p,black,
D 0 1p,black,
N 1
G 0.05i
''' % vars()
        # little padding for feature section
        h += 0.05

    # Custom text from .cfg file 
    if args.get('overlay_legend_text'):
        str = args.get('overlay_legend_text')
        txt += '''# custom text
L 9 Helvetica L %(str)s
G 0.05i
D 0 1p,black,
''' % vars()
        h += 0.25

    # mape scale 
    if args.get('overlay_legend_scale'): 
        txt += '''# map scale
G 0.05i
M 0 0 10000:km:l f -R%(R)s -J%(J)s
G 0.05i
''' % vars()
        h += 0.50

    # check for empty txt
    if txt == '':
        return

    # set args
    file_name = '%(ps)s.%(figure_id)s.pslegend.txt' % vars()
    file = open(file_name,'w')
    file.write('%s' % txt)
    file.close()
    args['pslegend_file'] = file_name

    dy = h + 0.2
    args['pslegend_D'] = '0/%(dy)f/5/%(h)f/BL' % vars()

    return file_name

#=====================================================================
# NOTE: sample strings for pslegend.txt files :
# S 0.1i v 0.25i/0.02/0.06/0.05 255/0/255 0.25p 0.3i This is a vector
## len/aw/hl/hw,
## arrowwidth/headlength/headwidth [Default is 0.075c/0.3c/0.25c (or 0.03i/0.12i/0.1i)]
#=====================================================================

#=====================================================================
def simple_plot_image( args ):
    ''' '''
    
#=====================================================================
def plot_image( args ):
    '''Simple sequence of commands to plot a grid file.'''

    # save any original -R passed in
    R = ''
    if args.has_key('R'):
        R = args['R']

    # start up the gmt seqence 
    start_page(args)

    # fill args with info on grid file; -R 
    grdinfo(args)

    # make cpt file from args info
    args['I'] = '' # invert default polar scale
    makecpt(args)

    # reset R back to what was passed in
    args['R'] = R 

    # keep the plot open 
    args['K'] = ' '
    args['O'] = ' '

    # plot the grid
    grdimage(args)

    # check for xy file 
    if args.get('xy'): 
        psxy(args)

    # close the page
    end_page(args)

    # clean up tmp files
    if args.get('C'):
        cpt = args['C']
        if cpt.startswith('tmp.') :
            cmd = 'rm -rf %(cpt)s' % vars()
            if verbose: print now(), 'plot_image: cmd =\n', cmd
            os.system(cmd)
#=====================================================================
#=====================================================================
def get_high_precision_region( dict ):
    '''call minmax on an xy or xyz file to set the region to a high
precision.  Helps when gridding xyz files by known increment.

Required arguments:
    dict['xyz_file]

Output arguments:
    dict['R']
    dict['grid_min']
    dict['grid_max'] 
'''
    xyz = dict['xyz_file']
    # query rotated xyz for -R values
    cmd = 'minmax --D_FORMAT='+'"%"'+'.12lg %(xyz)s' % vars() 
    pipe = os.popen(cmd)
    line = pipe.readline()
    Line=line.split('\t')
    Lon1=Line[1].strip('>')
    Lon2=Lon1.strip('<')
    Lat1=Line[2].strip('>')
    Lat2=Lat1.strip('<')
    # if minmax called on an xyz file also update grid_min
    # and gridmax
    R = Lon2+'/'+Lat2
    if len(Line) == 4: # condition NOT true for minmax of xy only
        dict['grid_min']=Line[3].strip('<').split('/')[0]
        dict['grid_max']=Line[3].strip('>\n').split('/')[1]
    else:
        R = R[:-2] # necessary correction for xy
    pipe.close()
    dict['R'] = R
    return R

#=====================================================================
#=====================================================================
def find_transition( dict ):
    '''Makes an xy file to plot a line at the transition zone
       in an x sect or y sect coor
       will make either transition_long.xy if x is set or
       transition_lat.xy if y is set.
       Transition zone depth is here hard wired at 670 km

Required in dictionary:
    1. dict['radius_inner']
    2. dict['z_lmantle']
    3. dict['nodez']
    4. dict['fi_min']
    5. dict['fi_max']
    6. dict['x']
    7. dict['y']
    8. dict['theta_max']
    9. dict['theta_min']


Resets the values in the dictionary:
    1. dict['xy']='transition_long.xy' or
       dict['xy']='transition_lat.xy'
    2. dict['N']='-N' (allows it to plot beyond the -R region.)

    '''
    r2d = 180.0/math.pi
    Transition_depth=float(670)
    z_mantle=float(dict['radius_inner'])
    z_transition=1-float(dict['z_lmantle'])
    prop=(1-z_transition)/(1-z_mantle)
    total_nodez=float(dict['nodez']) - 1 
    nodez_loc=int((total_nodez)-(prop*total_nodez) )
    print "================================="
    print "  nodez to 670 transition zone   "
    print    nodez_loc
    print "================================="
    # Transition zone along Latitudinal slice
    if dict['x'] <> 'none':
        longs=range(  int((float(dict['fi_min'])) *r2d) , int((float(dict['fi_max']))*r2d))
        print "transition zone at depth in km"
        print     Transition_depth
        print "transition zone at non-dimensional depth"
        print     z_transition 
        Trans_Long=open('transition_long.xy',"w")
        for i in range(len(longs)):
            Long=int(longs[i])	
            Trans_Long.write('%s %s \n' % (float(Long), z_transition))
        Trans_Long.close()
        dict['xy']='transition_long.xy'
        dict['N']='-N'
	# Transition zone along Longitudinal slice		
    if dict['y'] <> 'none':
        lats=range(  90 - int((float(dict['theta_max']))*r2d) , 90 - int((float(dict['theta_min'])) *r2d))
        print "transition zone at depth in km"
        print    Transition_depth
        print "transition zone at non-dimensional depth"
        print    z_transition 
        Trans_Lat=open('transition_lat.xy',"w")
        for i in range(len(lats)):
            Lat=int(lats[i])	
            Trans_Lat.write('%s %s\n' % (float(Lat), z_transition))
        Trans_Lat.close()
        dict['xy']='transition_lat.xy'
        dict['N']='-N'
	return

#=====================================================================
#=====================================================================
def random_velocities( dict ):
    '''Get a random selection of velocity vectors to plot

Required arguments:
 dict['velocity_file'] = name of velocity file produced by to create 
    
Optional arguments:
 dict['vel_density'] = default is 0.01

Output arguments
 dict['xy'] = 'sampled_velfile.xy'

Return value:
 none
'''
    # Define vel input file
    vel_infile=dict['velocity_file']
    print "============================="
    print "   Set Velocity file         "
    print "============================="
    print       vel_infile
    # Define output vel file
    print "============================="
    print "   Set Velocity out file     "
    print "============================="
    vel_outfile=open('sampled_velfile.xy','w')
    print       vel_outfile
    # parse velocity increment 
    vel_density = .01
    # get the sampling density from vel_density
    if dict.get('vel_density'): 
        vel_density = float(dict.get('vel_density'))
    print "============================="
    print " Using a sampling density of "
    print "============================="
    print   vel_density
    # How many entries are there in the vel file
    file_info=commands.getoutput("wc -l "+vel_infile+"")
    vel_file_length=int(file_info.split(' ')[0])
    sampling=int(vel_file_length*vel_density)
    # Store random numbers in a vector
    random_numbers=[]
    for i in range(sampling):
        random_numbers.append(random.randrange(1,vel_file_length,1))
    # Select vel vectors at random locations
    infile=open(vel_infile,'r')
    for i in range(vel_file_length):
	line=infile.readline()
	if random_numbers.count(i):
		vel_outfile.write(line)
    vel_outfile.close()
    infile.close()
    dict['velocity_sampled_file']='sampled_velfile.xy'
    dict['xy']='sampled_velfile.xy'
    return


#=====================================================================
#=====================================================================
def xyz2grd_smoothing( dict ) :
    '''Create a grid file using GMT calls. Same as xyz2grd, however 
    this one filters the grid according to some filter size. This is
    particularly useful for z-sections where sometimes surface results
    in undefined values. So in addition to the requirements for xyz2grd
    user also needs to specify the value ['width'] in the dictionary.
    Uses a gaussian filter -Fg with -D2 and the gmt command grdfilter.
Requires:
    1. dict['xyz_file'] input data
    2. dict['grid_increment'] sets -I value
    3. dict['grid_min'] sets -Ll value
    4. dict['grid_max'] sets -Lu value
    5. dict['grid_tension'] sets -T value
    6. dict['headers'] sets -H value
    5. dict['file_in_latlon'] sets -: flag
    5. dict['width'] sets -T value
Returns
    1. grid_file

'''
    # file parameters
    xyz_file = dict['xyz_file']
    mean_file="%s_median.xyz" % xyz_file
    grid_file = "%s.grd" % xyz_file
    grid_file2 = "%s.grd" % xyz_file
    # geography
    R = dict['R']
    # filter width
    width = float(dict['width'])
    # grid parameters
    grid_increment = dict['grid_increment']
    grid_min = dict['grid_min']
    grid_max = dict['grid_max']
    grid_tension = dict['grid_tension']
    # smooth the data into a mean file 
    cmd = 'blockmedian %(xyz_file)s ' % vars()
    if dict.get('file_headers'): 
        H = dict.get('file_headers')
        cmd += '-H%(H)s ' % vars()
    if dict.get('file_in_latlon'): cmd += '-: ' % vars()
    cmd += '-I%(grid_increment)s -R%(R)s > %(mean_file)s' % vars()
    os.system(cmd)
    # check for file 
    if not os.path.exists( mean_file ):
        msg = 'blockmean may have failed: file not found: %(mean_file)s' % vars()
        raise ValueError, msg
    # create the grid
    cmd = 'surface %(mean_file)s ' % vars()
    if dict.get('file_headers'): 
        H = dict.get('file_headers')
        cmd += '-H%(H)s ' % vars()
    if dict.get('file_in_latlon'): cmd += '-: ' % vars()
    if grid_min != 'none' or grid_max != 'none':
        cmd += '-Ll%(grid_min)s -Lu%(grid_max)s ' % vars()
    cmd += '-G%(grid_file2)s -I%(grid_increment)s -R%(R)s -T%(grid_tension)s' % vars()
    os.system(cmd)
    # check for file 
    if not os.path.exists( grid_file ):
        msg = 'surface may have failed: file not found: %(mean_file)s' % vars()
        raise ValueError, msg
    # filter the grid file
    cmd = 'grdfilter %(grid_file2)s -D2 -Fg%(width)g -V -G%(grid_file)s' % vars()
    os.system(cmd)   
    # clean up
    cmd = "rm -rf %(mean_file)s" % vars()
    os.system(cmd)
    return grid_file

#=====================================================================
#=====================================================================

def test_plot( argv ):
    '''test module functions'''

    global verbose
    verbose = True
    # establish empty dict
    args = {}

    # fill dict with cmd line data
    grid = argv[1]
    args['grid'] = grid

    # fill dict with test data 
    args['ps'] = '_test_plot_%(grid)s.ps' % vars()

    # fill dict with user choices
    #args['page_orientation'] = 'landscape'
    #args['J'] = 'H0/9'
    #args['X'] = 'c'
    #args['Y'] = 'c'

    args['page_orientation'] = 'portrait'
    args['B'] = 'a30f30:.%(grid)s:WESN' % vars()
    args['J'] = 'H0/7'
    args['J'] = 'M7'
    args['X'] = 'c'
    args['Y'] = 'c'
   
    args['dz'] = '1' # grid z units

    args['page_title'] = 'test title string'
    args['page_number'] = 'test page number string'

    # start up the gmt seqence 
    start_page(args)

    # fill args with info on grid file; -R 
    grdinfo(args)

    # make cpt file from args info
    args['I'] = '' # invert default polar scale
    makecpt(args)

    # keep the plot open 
    args['K'] = ' '
    args['O'] = ' '

    # plot the grid
    grdimage(args)

    Core_Util.tree_print(args)

    end_page(args)
#=====================================================================
#=====================================================================
def makecpt_def( argv ):
    '''This function prints to std out data suitable for a .cpt file to highlight deformation color values'''

    # Going fwd in time, triangles that are expanding will have positive dilitation 
    # and a negative color value, for example:
    # colour_num = -14.8281 ; dilitation = 1.48547e-15
    # color these triangles shades of red 

    # Going fwd in time, triangles that are shrinking will have negative dilitation
    # and a positive color value: for example:
    # colour_num = 14.9088 ; dilitation = -1.23372e-15
    # color these triangles shades of blue 

    # set the header lines
    header = '''#
# cpt file for deformation tests
#'''

    # arrays to hold table entries for red side and blue side
    r_lines = []
    b_lines = []

    #
    # create the red side
    #
    r_lower = -20
    r_upper = -13
    r_inc = .1

    # total number of entries in table
    r_count  = ( r_upper - r_lower) / r_inc

    # color values range from 0 to 255 
    # with increment dependant on number of entries 
    c_lower = 0
    c_upper = 255
    c_inc = 255 / r_count 

    # loop over the index values and build the red lines
    i = 0
    while i < r_count:
        vi = r_lower + (i * r_inc );
        vn = vi + r_inc;

        # as data values range from r_lower to r_upper,
        # colors range from red to white 
        r = c_upper
        b = c_upper - (i * c_inc) 
        g = c_upper - (i * c_inc)

        line = "%(vi)g\t%(r)d\t%(g)d\t%(b)d\t%(vn)g\t%(r)d\t %(g)d\t%(b)d" % vars()
        r_lines.append(line)

        # update the loop
        i = i + 1;
    # end of loop 

    # add a line for values from r_upper to 0.0; color them white 
    r = g = b = 255;
    line = "%(r_upper)g\t%(r)d\t%(g)d\t%(b)d\t0.0\t%(r)d\t %(g)d\t%(b)d" % vars()
    r_lines.append(line)

    # 
    # Create Blue side
    # 
    b_lower = 13
    b_upper = 20
    b_inc = .1

    # total number of entries on blue side of table
    b_count  = ( b_upper - b_lower) / b_inc

    # color values range from 0 to 255 
    # with increment dependant on number of entries 
    c_lower = 0
    c_upper = 255
    c_inc = 255 / b_count 

    # Add a line to account for values from 0 to b_lower 
    r = g = b = c_upper
    line = "0.0\t%(r)d\t%(g)d\t%(b)d\t%(b_lower)g\t%(r)d\t %(g)d\t%(b)d" % vars()
    b_lines.append(line)

    # loop over the index values and build the lines
    i = 0
    while i < b_count:
        vi = b_lower + (i * b_inc );
        vn = vi + b_inc;

        # as data values range from b_lower to b_upper,
        # colors range from blue to white
        r = c_lower + (i * c_inc)
        g = c_lower + (i * c_inc)
        b = c_upper 

        line = "%(vi)g\t%(r)d\t%(g)d\t%(b)d\t%(vn)g\t%(r)d\t %(g)d\t%(b)d" % vars()
        b_lines.append(line)

        # update the loop
        i = i + 1;
    # end of loop

    # any values outside the ranges of r_lower to r_upper and b_lower to b_upper will be white 
    footer = '''B\t255\t255\t255
F\t255\t255\t255
N\t255\t255\t255'''

    print header
    for l in r_lines:
        print l
    for l in b_lines:
        print l
    print footer

#=====================================================================
#=====================================================================
if __name__ == "__main__":
    import Core_GMT
    if len(sys.argv) > 1:
        # test_plot( sys.argv )
        makecpt_def( sys.argv );
    else:
      help(Core_GMT)
#=====================================================================
# to create a deformation color scale run this command:
#
# $ Core_GMT.py 1 > /private/tmp/deformation_test.cpt
#
#=====================================================================
#=====================================================================

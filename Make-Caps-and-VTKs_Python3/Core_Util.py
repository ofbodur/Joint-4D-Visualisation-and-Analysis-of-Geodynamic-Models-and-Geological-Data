#!/usr/bin/python
#=====================================================================
#                      Python Scripts for CitComS 
#         Preprocessing, Data Assimilation, and Postprocessing
#                  ---------------------------------
#             (c) California Institute of Technology 2008
#                        ALL RIGHTS RESERVED
#=====================================================================
''' A set of general purpose functions.

This module has functions to read and parse various file types associated 
with CitComS, GPlates, the plot_page system, and other data sets.

Most functions take a dictionary of arguments to pass and adjust paramters.
'''
#=====================================================================
#=====================================================================
# imports
import getopt, os, string, sys, math, time, datetime, pprint
#=====================================================================
# Global variables
verbose = False
#=====================================================================
#=====================================================================
def now():
    '''redefine now() with short format yyyy-mm-dd hh:mm:ss'''
    return str(datetime.datetime.now())[:19]
#=====================================================================
def tree_print(arg):
    '''print arg in tree form'''
    pprint.PrettyPrinter(indent=2).pprint(arg)
#=====================================================================
def command(command):
    '''call a unix command and return the stdout as a string'''
    child = os.popen(command)
    data = child.read()
    err = child.close()
    if err:
        raise (RuntimeError, '%s failed w/ exit code %d') % (command, err)
    return data
#=====================================================================
#=====================================================================
def test( argv ):
    '''self test'''
    global verbose
    verbose = True 

    file = argv[1]
    print (now(), 'test: file = %(file)s') % vars()
#=====================================================================
#=====================================================================
if __name__ == "__main__":
    import Core_Util

    if len( sys.argv ) > 1:
        # process sys.arv as file names for testing 
        test( sys.argv )
    else:
        # print module documentaiton and exit
        help(Core_Util)
#=====================================================================
#=====================================================================
def sample_function( args ):
    '''sample_function shows how to code a simple function.

Required arguments:
 x = some number
 y = some other number
    
Optional arguments:
 none

Output arguments:
 z = set to x + y

Return value:
 w = set to x * y
'''

    # parse args dictionary
    x = args['x']
    y = args['y']

    # process data
    z = x + y
    w = x * y

    # fill args dictionary 
    args['z'] = z

    # return value
    return w
#=====================================================================
#=====================================================================
def test_sample( argv ):
    '''test calls to sample_function()'''
    global verbose
    verbose = True 

    # working dictionary 
    dict = {} 

    #
    # test with fixed input
    #
    x = 3
    y = 4
    dict['x'] = x
    dict['y'] = y
    print (now(), 'test: x = %(x)s') % vars()
    print (now(), 'test: y = %(y)s') % vars()

    # call function, get new values 
    w = sample_function (dict) 
    z = dict['z'] 
    print (now(), 'test: z = %(z)s') % vars()
    print (now(), 'test: w = %(w)s') % vars()

    # raise and error message if test fails:
    if w != x * y: 
       msg = 'w is not correct, for x=', x, 'and y=', y
       raise (ValueError, msg)

    if z != x + y:
       msg = 'w is not correct, for x=', x, 'and y=', y
       raise (ValueError, msg)


    # 
    # test with user input input, 
    #
    x = int ( argv[1] ) 
    y = int ( argv[2] )
    dict['x'] = x
    dict['y'] = y
    print (now(), 'test: x = %(x)s') % vars()
    print (now(), 'test: y = %(y)s') % vars()

    # call function, get new values 
    w = sample_function (dict) 
    z = dict['z'] 
    print (now(), 'test: z = %(z)s') % vars()
    print (now(), 'test: w = %(w)s') % vars()

    # raise and error message if test fails:
    if w != x * y:
       msg = 'w is not correct, for x=', x, 'and y=', y
       raise (ValueError, msg)

    if z != x + y:
       msg = 'w is not correct, for x=', x, 'and y=', y
       raise (ValueError, msg)
#=====================================================================
#=====================================================================

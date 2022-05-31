#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#<LicenseText>
#
# CitcomS.py by Eh Tan, Eun-seo Choi, and Pururav Thoutireddy.
# Copyright (C) 2002-2005, California Institute of Technology.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#</LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

'''Automatically find the input parameters from CitcomS input file and run
'batchcombine.py'

usage: autocombine.py machinefile inputfile step1 [step2 [...] ]
'''

# default values for CitcomS input
defaults = {'output_format': 'ascii',
            'output_optional': 'surf,botm,pressure,stress,heating,composition',
            'buoy_type': 1,
            'tracer_flavors': 2,
            'nprocx': 1,
            'nprocy': 1,
            'nprocz': 1,
            'nodex': 9,
            'nodey': 9,
            'nodez': 9}



def normalize_optional(output_optional):
    #fields = ['surf','botm','horiz_avg','pressure','stress']
    fields=[]
    for opt in output_optional.split(','):
        ## remove the leading/trailing whitespaces
        opt = opt.strip()

        ## retain fields that are node-based
        if opt in ('pressure','stress','temperature'):
            fields.append(opt)


    return ','.join(fields)



if __name__ == '__main__':

    import sys
    import batchcombine as bc

    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)

    machinefile = sys.argv[1]
    inputfile = sys.argv[2]
    timesteps = [int(i) for i in sys.argv[3:]]

    # parse input
    from parser import Parser
    parser = Parser(defaults)
    parser.read(inputfile)

    datadir = parser.getstr('datadir')
    datafile = parser.getstr('datafile')
    output_format = parser.getstr('output_format')

    if output_format != 'ascii':
        print("Error: don't know how to combine the output", \
              "(output_format=%s)" % output_format)
        sys.exit(1)

    output_optional = parser.getstr('output_optional')
    optional_fields = normalize_optional(output_optional)

    buoy_type = parser.getint('buoy_type')
    nflavors = parser.getint('tracer_flavors')
    if buoy_type == 0:
        ncompositions = nflavors
    elif buoy_type == 1:
        ncompositions = nflavors - 1
    else:
        print("Error: don't know how to combine the output", \
              "(buoy_type=%d)" % buoy_type)
        sys.exit(1)



    nodex = parser.getint('nodex')
    nodey = parser.getint('nodey')
    nodez = parser.getint('nodez')
    ncap = parser.getint('nproc_surf')
    nprocx = parser.getint('nprocx')
    nprocy = parser.getint('nprocy')
    nprocz = parser.getint('nprocz')

    totalnodes = nprocx * nprocy * nprocz * ncap
    nodelist = bc.machinefile2nodelist(machinefile, totalnodes)

    for timestep in timesteps:
        # combining coord, velo, temp, and visc
        bc.batchcombine(nodelist, datadir, datafile, timestep,
                        nodex, nodey, nodez,
                        ncap, nprocx, nprocy, nprocz,
                        'coord,velo,visc,comp_nd',0)
                        #'coord,velo,visc,pressure,stress,comp_nd', 0)
        print(optional_fields)
        # combining optional fields, if necessary
        #if optional_fields:
        #    bc.batchcombine(nodelist, datadir, datafile, timestep,
        #                    nodex, nodey, nodez,
        #                    ncap, nprocx, nprocy, nprocz,
        #                    optional_fields, ncompositions)



# version
# $Id: autocombine.py 7840 2007-08-17 18:33:36Z tan2 $

# End of file

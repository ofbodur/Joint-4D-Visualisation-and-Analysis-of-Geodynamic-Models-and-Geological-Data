#!/usr/bin/env python2.7 
# Update: Converted to python 3.8 on 30 May 2022 (Ömer F. Bodur), University of Wollongong.

'''
/*
 * =====================================================================================
 * Copyright (C) 2014 Rakib Hassan (rakib.hassan@sydney.edu.au)
 *
 * This program is free software; you can redistribute it and/or modify it under 
 * the terms of the GNU General Public License as published by the Free Software 
 * Foundation; either version 2 of the License, or (at your option) any later 
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple 
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * ===================================================================================== 

 * ===================================================================================== 
 * Updated by Ömer F. Bodur at University of Wollongong (30/05/2022).

 * This program is free software; you can redistribute it and/or modify it under 
 * the terms of the GNU General Public License as published by the Free Software 
 * Foundation; either version 3 of the License, or (at your option) any later 
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple 
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * ===================================================================================== 
 */
'''

import vtk
import os
import sys
import math
import getopt, os, string, sys, math, time, subprocess
import datetime
import pprint

import CitcomParser, Core_Util, Core_Citcom
from scipy import spatial
import numpy

msg=\
"""

Extracts data from a regional cap file, converts coords to cartesian and 
outputs a .vtu file.
Usage: ./capToVtk <capfile> <pidfile>
"""



# import h5py
# from vtk.util.vtkAlgorithm import VTKPythonAlgorithmBase
# from vtk.numpy_interface import dataset_adapter as dsa
 
# class HDF5Writer(VTKPythonAlgorithmBase):
#     def __init__(self):
#         VTKPythonAlgorithmBase.__init__(self,
#             nInputPorts=1, inputType='vtkUnstructuredGrid',
#             nOutputPorts=0)
#
#         self.__FileName = ""
#         self.__NumberOfPieces = 1
#         self.__CurrentPiece = 0
#
#     def RequestData(self, request, inInfo, outInfo):
#         info = inInfo[0].GetInformationObject(0)
#         inp = dsa.WrapDataObject(vtk.vtkDataSet.GetData(info))
#
#         if self.__CurrentPiece == 0:
#               self.__File = h5py.File(self.__FileName, 'w')
#
#         grp = self.__File.create_group("piece%d" % self.__CurrentPiece)
#         grp.attrs['bounds'] = inp.GetBounds()
#
#         grp.create_dataset("cells", data=inp.Cells)
#         grp.create_dataset("cell_types", data=inp.CellTypes)
#         grp.create_dataset("cell_locations", data=inp.CellLocations)
#
#         grp.create_dataset("points", data=inp.Points)
#
#         pdata = grp.create_group("point_data")
#         for name in inp.PointData.keys():
#             pdata.create_dataset(name, data=inp.PointData[name])
#
#         if self.__CurrentPiece < self.__NumberOfPieces - 1:
#             # If we are not done, ask the pipeline to re-execute us.
#             self.__CurrentPiece += 1
#             request.Set(
#                 vtk.vtkStreamingDemandDrivenPipeline.CONTINUE_EXECUTING(),
#                 1)
#         else:
#             # Stop execution
#             request.Remove(
#                 vtk.vtkStreamingDemandDrivenPipeline.CONTINUE_EXECUTING())
#             self.__File.close()
#             del self.__File
#         return 1
#
#     def RequestInformation(self, request, inInfo, outInfo):
#         # Reset values.
#         self.__CurrentPiece = 0
#         return 1
#
#     def RequestUpdateExtent(self, request, inInfo, outInfo):
#         info = inInfo[0].GetInformationObject(0)
#         info.Set(
#             vtk.vtkStreamingDemandDrivenPipeline.UPDATE_NUMBER_OF_PIECES(),
#             self.__NumberOfPieces)
#         info.Set(
#             vtk.vtkStreamingDemandDrivenPipeline.UPDATE_PIECE_NUMBER(),
#             self.__CurrentPiece)
#         return 1
#
#     def SetFileName(self, fname):
#         if fname != self.__FileName:
#             self.Modified()
#             self.__FileName = fname
#
#     def GetFileName(self):
#         return self.__FileName
#
#     def SetNumberOfPieces(self, npieces):
#         if npieces != self.__NumberOfPieces:
#             self.Modified()
#             self.__NumberOfPieces = npieces
#
#     def GetNumberOfPieces(self):
#         return self.__NumberOfPieces




# # Dimensionalize variables
# # Model parameters
# SurfTemp=0.08
# alpha0= 3.0e-5# thermal expansivity
# rho0= 4000.0 # ref. density Nico had 3340
# g0= 9.81 #grav. acc.
# # DeltaT= #Temperature diff.
# k0= 1.0e-6 #thermal diff. m^2/s
radiusEarth= 6371.0e3 # radius of the Earth (metres)
# eta0= 1.0e21 #ref. viscosity
# Ra=8.6e08
#
# PressureDim=(eta0*k0)/(radiusEarth**2)
# # Ra=(alpha0*rho0*g0*DeltaT*((hm)**3))/(k0*eta0)
#
# DeltaT=(Ra*k0*eta0)/(alpha0*rho0*g0*(radiusEarth**3))
# print(DeltaT*SurfTemp -273., ' is the surface temperature in C')
#
# # TempDim=DeltaT*(1+SurfTemp)
#
# velocityDimFactor=(3.154e+9)*k0/radiusEarth # cm/yr


def rtp2xyz(r, theta, phi):
    rst = r * math.sin(theta)
    xout = [0]*3
    xout[0] = rst * math.cos(phi)	    # x 
    xout[1] = rst * math.sin(phi) 	    # y 
    xout[2] = r * math.cos(theta)       # z

    return xout

    
def vrtp2vxyz(vr, vtheta, vphi, r, theta, phi):
    sinTheta = math.sin(theta)
    cosTheta = math.cos(theta)
    
    sinPhi = math.sin(phi)
    cosPhi = math.cos(phi)

    xout=[0]*3

    xout[0] = (sinTheta*cosPhi*vr + cosTheta*cosPhi*vtheta - sinPhi*vphi)
    xout[1] = (sinTheta*sinPhi*vr + cosTheta*sinPhi*vtheta + cosPhi*vphi)
    xout[2] = (cosTheta*vr - sinTheta*vtheta)

    return xout  
    
def combineCompos(compos1,compos2,compos3,compos4):
    xout=[0]*4
    xout[0]=compos1
    xout[1]=compos2
    xout[2]=compos3
    xout[3]=compos4
    
    return xout
    
  
def processData(capFile, paramsDict):
    nodez = int(paramsDict['nodez'])
    nodey = int(paramsDict['nodey'])
    nodex = int(paramsDict['nodex'])

    pnodes=[]
    # pvelocity=[]
    pNDVelocity=[]
    ptemperature=[]
    pviscossity=[]
    # pstress=[]
    # ppressure=[]
    # pcomposition=[]
    pvr=[]
    pradius=[]
    pdepth=[]
    pcompos=[]
    pLLSVP=[]
    #read data
    i=0
    f = open(capFile)
    for line in f:
        i = i + 1 #skip header
        if (i==1): continue
        
        cols = line.split(" ")

        theta   = float(cols[0]) #colatitude
        phi     = float(cols[1]) #longitude
        r       = float(cols[2]) #radius
        vtheta  = float(cols[3])
        vphi    = float(cols[4])
        vr      = float(cols[5])
        temp    = float(cols[6])
        visc    = float(cols[7])
        compos_1 = float(cols[8])
        compos_2 = float(cols[9])
        compos_3 = float(cols[10])
        compos_4 = float(cols[11])
        compos_5 = float(cols[12])
        
        # stheta_theta = float(cols[8])
        #         sphi_phi = float(cols[9])
        #         sr_r = float(cols[10])
        #         stheta_phi = float(cols[11])
        #         Stheta_r = float(cols[12])
        #         Sphi_r = float(cols[13])
        #         pres = float(cols[14])
       
        
        xyz     = rtp2xyz(r, theta, phi)
        # vxyz    = vrtp2vxyzDimensional(vr, vtheta, vphi, r, theta, phi)
        vxyzND  = vrtp2vxyz(vr, vtheta, vphi, r, theta, phi)
        # sAsis   = takeStressAsItIs(stheta_theta, sphi_phi, sr_r, stheta_phi, Stheta_r, Sphi_r)
        pnodes.append(xyz)
        # pvelocity.append(vxyz)
        pNDVelocity.append(vxyzND)
        ptemperature.append(temp)
        pviscossity.append(visc)
        composConnected=combineCompos(compos_1,compos_2,compos_3,compos_4)
        
        pcompos.append(composConnected)
        pLLSVP.append(compos_5)
        
        # pstress.append(sAsis)
        #         ppressure.append(pres)
        #         pcomposition.append(compos)
        pvr.append(vr) # radial vel.
        pradius.append(r)
        pdepth.append(1.0-r)
        
         
    #end for
    phoriz_avg=[]
    pradius_horiz_avg=pradius[0:65]
    pradialheatadv=[]
    pVrms_horiz=[]
    pVrms_vertical=[]
    pComp_avg_1=[]
    pComp_avg_2=[]
    pComp_avg_3=[]
    pComp_avg_4=[]
    pComp_avg_5=[]    
    print("burada once")
    # print pdepth
    i=0
    # if len(capfileFile)==71:
    tstep=capFile[47:]
    print(capFile)
    print(tstep)    
   
    f2 = open('/Volumes/Accelsior4M2/gld422/HorizontalAverage_opt/gld422.horiz_avg.0.'+tstep)
    
    f3 = open('/Volumes/Accelsior4M2/gld422/HorizontalAverage_opt/gld422.horiz_avg.1.'+tstep)
    
    for line in f2:
        # if (i==0): continue
        columns = line.split(" ")
        # print line
        # print "burada"
        
        # radius_horiz_avg  = float(columns[0])
        horiz_avg = float(columns[1])
        Vrms_horizAdd= float(columns[2])
        Vrms_verticalAdd = float(columns[3])
        Comp_avg_1= float(columns[4])
        Comp_avg_2= float(columns[5])
        Comp_avg_3= float(columns[6])
        Comp_avg_4= float(columns[7])
        Comp_avg_5= float(columns[8])
        
        # pradius_horiz_avg.append(radius_horiz_avg)
        phoriz_avg.append(horiz_avg)
        pVrms_horiz.append(Vrms_horizAdd)
        pVrms_vertical.append(Vrms_verticalAdd)
        pComp_avg_1.append(Comp_avg_1)
        pComp_avg_2.append(Comp_avg_2)
        pComp_avg_3.append(Comp_avg_3)
        pComp_avg_4.append(Comp_avg_4)
        pComp_avg_5.append(Comp_avg_5)
    ll=0            
    for line in f3:
        columns = line.split(" ")
        # print line
        # print "burada"
        
        # radius_horiz_avg  = float(columns[0])
        horiz_avg = float(columns[1])
        Vrms_horizAdd= float(columns[2])
        Vrms_verticalAdd = float(columns[3])
        Comp_avg_1= float(columns[4])
        Comp_avg_2= float(columns[5])
        Comp_avg_3= float(columns[6])
        Comp_avg_4= float(columns[7])
        Comp_avg_5= float(columns[8])
        
        # pradius_horiz_avg.append(radius_horiz_avg)
        if ll>0:
            phoriz_avg.append(horiz_avg)  
            pVrms_horiz.append(Vrms_horizAdd)
            pVrms_vertical.append(Vrms_verticalAdd)
            pComp_avg_1.append(Comp_avg_1)
            pComp_avg_2.append(Comp_avg_2)
            pComp_avg_3.append(Comp_avg_3)
            pComp_avg_4.append(Comp_avg_4)
            pComp_avg_5.append(Comp_avg_5)
        ll=ll+1
        
 
        
    # print phoriz_avg
    # print(phoriz_avg)
    tree = spatial.cKDTree(pnodes, 16)

    connectivity = []    
    i=1
    j=1
    k=1
    for n in range(nodex*nodey*nodez - (nodex*nodez)):
        if ((i%nodez)==0):   #X-Values
            j=j+1             #Count Y-Values
        #end if
        
        if ((j%nodex)==0):
            k=k+1             #Count Z-Values
        #end if
							
        if (((i%nodez) != 0) and ((j%nodex) != 0)):            #Check if Box can be created
            #Create Connectivity
            cell = [0]*8
            cell[0] = n                                #1
            cell[1] = (cell[0]+nodez)                  #2
            cell[2] = (cell[1]+nodez*nodex)            #3
            cell[3] = (cell[0]+nodez*nodex)            #4
            cell[4] = (cell[0]+1)                      #5
            cell[5] = (cell[4]+nodez)                  #6
            cell[6] = (cell[5]+nodez*nodex)            #7
            cell[7] = (cell[4]+nodez*nodex)            #8
            
            connectivity.append(cell)
        #end if
        
        i=i+1
    #end for
    
    math = vtk.vtkMath()
    points = vtk.vtkPoints()
    for i in range(len(pnodes)):
        points.InsertNextPoint(pnodes[i])
    #end for

    npoints = points.GetNumberOfPoints()
    
    NonDimVelocity    = vtk.vtkFloatArray()
    NonDimVelocity.SetName('NonDimVelocity')
    NonDimVelocity.SetNumberOfComponents(3)
    NonDimVelocity.SetNumberOfTuples(npoints)


    NonDimTemperature = vtk.vtkFloatArray()
    NonDimTemperature.SetName('NonDimTemperature')
    NonDimTemperature.SetNumberOfValues(npoints)
    
 
    NonDimViscosity   = vtk.vtkFloatArray()
    NonDimViscosity.SetName('NonDimViscosity')
    NonDimViscosity.SetNumberOfValues(npoints)
    
    Composition   = vtk.vtkFloatArray()
    Composition.SetName('Composition')
    Composition.SetNumberOfComponents(4)
    Composition.SetNumberOfValues(npoints)
    
    LLSVp = vtk.vtkFloatArray()
    LLSVp.SetName('LLSVP')
    LLSVp.SetNumberOfValues(npoints)
    
    
    # stress    = vtk.vtkFloatArray()
    #     stress.SetName('Stress')
    #     stress.SetNumberOfComponents(6)
    #     stress.SetNumberOfTuples(npoints)
    #
    #     pressure   = vtk.vtkFloatArray()
    #     pressure.SetName('Pressure(Pa)')
    #     pressure.SetNumberOfValues(npoints)
    #
    #     composition   = vtk.vtkFloatArray()
    #     composition.SetName('Composition')
    #     composition.SetNumberOfValues(npoints)
    #
    #     radialvelocity   = vtk.vtkFloatArray()
    #     radialvelocity.SetName('Vr(cm/yr)')
    #   radialvelocity.SetNumberOfValues(npoints)
    #  
    NonDimRadialVelocity = vtk.vtkFloatArray()
    NonDimRadialVelocity.SetName('NonDimRadialVelocity')
    NonDimRadialVelocity.SetNumberOfValues(npoints)
    
    
    radius   = vtk.vtkFloatArray()
    radius.SetName('Radius(km)')
    radius.SetNumberOfValues(npoints)
    
    NonDimRadius = vtk.vtkFloatArray()
    NonDimRadius.SetName('NonDimRadius')
    NonDimRadius.SetNumberOfValues(npoints)
    
    depth = vtk.vtkFloatArray()
    depth.SetName('Depth(km)')
    depth.SetNumberOfValues(npoints)    
    
    NonDimDepth = vtk.vtkFloatArray()
    NonDimDepth.SetName('NonDimDepth')
    NonDimDepth.SetNumberOfValues(npoints)
    

    # temperature_Anomaly= vtk.vtkFloatArray()
 #    temperature_Anomaly.SetName('Temperature_Anomaly(K)')
 #    temperature_Anomaly.SetNumberOfValues(npoints)
 # 
    NonDimTemperature_Anomaly= vtk.vtkFloatArray()
    NonDimTemperature_Anomaly.SetName('NonDimTemperature_Anomaly')
    NonDimTemperature_Anomaly.SetNumberOfValues(npoints)
    
    # radialheatadv= vtk.vtkFloatArray()
  #   radialheatadv.SetName('Radial_Heat_Advection(mK/yr)')
  #   radialheatadv.SetNumberOfValues(npoints)
  #  
    # horizontal_av= vtk.vtkFloatArray()
 #    horizontal_av.SetName('Horiz_Avg_Temp(K)')
 #    horizontal_av.SetNumberOfValues(npoints)
 # 
    NonDimHorizontal_Av=vtk.vtkFloatArray()
    NonDimHorizontal_Av.SetName('NonDimHorizontal_Av')
    NonDimHorizontal_Av.SetNumberOfValues(npoints)
    
    NonDim_Vrms_horiz = vtk.vtkFloatArray()
    NonDim_Vrms_horiz.SetName('NonDim_Vrms_horiz')
    NonDim_Vrms_horiz.SetNumberOfValues(npoints)
    
    NonDim_Vrms_radial = vtk.vtkFloatArray()
    NonDim_Vrms_radial.SetName('NonDim_Vrms_radial')
    NonDim_Vrms_radial.SetNumberOfValues(npoints)
    
    Compositional_Av1=vtk.vtkFloatArray()
    Compositional_Av1.SetName('Comp_avg_1')
    Compositional_Av1.SetNumberOfValues(npoints)
    
    Compositional_Av2=vtk.vtkFloatArray()
    Compositional_Av2.SetName('Comp_avg_2')
    Compositional_Av2.SetNumberOfValues(npoints)
    
    Compositional_Av3=vtk.vtkFloatArray()
    Compositional_Av3.SetName('Comp_avg_3')
    Compositional_Av3.SetNumberOfValues(npoints)
    
    Compositional_Av4=vtk.vtkFloatArray()
    Compositional_Av4.SetName('Comp_avg_4')
    Compositional_Av4.SetNumberOfValues(npoints)
    
    Compositional_Av5=vtk.vtkFloatArray()
    Compositional_Av5.SetName('Comp_avg_5')
    Compositional_Av5.SetNumberOfValues(npoints)
    
    
    
    # TempAnomalies_5Cases=vtk.vtkFloatArray()
    #     TempAnomalies_5Cases.SetName('TempAnomalies_5Cases')
    #     TempAnomalies_5Cases.SetNumberOfValues(npoints)
    
    Hot_Rising_Anomalies=vtk.vtkFloatArray()
    Hot_Rising_Anomalies.SetName('Hot_Rising_Anomalies')
    Hot_Rising_Anomalies.SetNumberOfValues(npoints)
    
    Hot_Sinking_Anomalies=vtk.vtkFloatArray()
    Hot_Sinking_Anomalies.SetName('Hot_Sinking_Anomalies')
    Hot_Sinking_Anomalies.SetNumberOfValues(npoints)
    
    Cold_Rising_Anomalies=vtk.vtkFloatArray()
    Cold_Rising_Anomalies.SetName('Cold_Rising_Anomalies')
    Cold_Rising_Anomalies.SetNumberOfValues(npoints)
    
    Cold_Sinking_Anomalies=vtk.vtkFloatArray()
    Cold_Sinking_Anomalies.SetName('Cold_Sinking_Anomalies')
    Cold_Sinking_Anomalies.SetNumberOfValues(npoints)
    # Calculator for shape file: 
    #cos(coordsY*3.14/180)*cos(coordsX*3.14/180)*iHat + cos(coordsY*3.14/180)*sin(coordsX*3.14/180)*jHat +  sin(coordsY*3.14/180)*kHat
    
    for i in range(npoints):
        Hot_Rising_Anomalies.SetValue(i,0.0)
        Hot_Sinking_Anomalies.SetValue(i,0.0)
        Cold_Rising_Anomalies.SetValue(i,0.0)
        Cold_Sinking_Anomalies.SetValue(i,0) 
    
    # print pradius_horiz_avg
    for i in range(npoints):
        # velocity.InsertTuple3(i, pvelocity[i][0], pvelocity[i][1], \
        #                                  pvelocity[i][2])
        #
        NonDimVelocity.InsertTuple3(i, pNDVelocity[i][0], pNDVelocity[i][1], \
                                 pNDVelocity[i][2])      
                                                                                   
        # temperature.SetValue(i, DeltaT*(ptemperature[i]+SurfTemp))
        
        NonDimTemperature.SetValue(i,ptemperature[i])
        
        
        # if( (ptemperature[i]-(phoriz_avg[pradius_horiz_avg.index(pradius[i])]))>0.0 and pvr[i]*velocityDimFactor>0.1 ):
        if ((ptemperature[i]-phoriz_avg[pradius_horiz_avg.index(pradius[i])])>0.0 and pvr[i]>0.0 ):    
            Hot_Rising_Anomalies.SetValue(i,(ptemperature[i]-phoriz_avg[pradius_horiz_avg.index(pradius[i])])*pvr[i])
        else:
            Hot_Rising_Anomalies.SetValue(i,0.0)
            
        
        # if( (ptemperature[i]-(phoriz_avg[pradius_horiz_avg.index(pradius[i])]))>0.0 and pvr[i]*velocityDimFactor<-0.1 ):
        if ((ptemperature[i]-phoriz_avg[pradius_horiz_avg.index(pradius[i])])>0.0 and pvr[i]<0.0 ):    
        
            Hot_Sinking_Anomalies.SetValue(i,(ptemperature[i]-phoriz_avg[pradius_horiz_avg.index(pradius[i])])*pvr[i]*(-1))
        else:
            Hot_Sinking_Anomalies.SetValue(i,0.0)
            
            
        # if( (ptemperature[i]-(phoriz_avg[pradius_horiz_avg.index(pradius[i])]))<0.0 and pvr[i]*velocityDimFactor>0.1 ):
        if ((ptemperature[i]-phoriz_avg[pradius_horiz_avg.index(pradius[i])])<0.0 and pvr[i]>0.0 ):    
        
            # TempAnomalies_5Cases.SetValue(i,3.0) # cold rising
            Cold_Rising_Anomalies.SetValue(i,(ptemperature[i]-phoriz_avg[pradius_horiz_avg.index(pradius[i])])*pvr[i]*(-1))
            
        else:
            Cold_Rising_Anomalies.SetValue(i,0.0)     
            
               
        # if( (ptemperature[i]-(phoriz_avg[pradius_horiz_avg.index(pradius[i])]))<0.0 and pvr[i]*velocityDimFactor<-0.1 ):
        if ((ptemperature[i]-phoriz_avg[pradius_horiz_avg.index(pradius[i])])<0.0 and pvr[i]<0.0 ):    
        
            # TempAnomalies_5Cases.SetValue(i,4.0) # cold sinking
            Cold_Sinking_Anomalies.SetValue(i,(ptemperature[i]-phoriz_avg[pradius_horiz_avg.index(pradius[i])])*pvr[i])
        else:
            # TempAnomalies_5Cases.SetValue(i,0.0) # stationary
            Cold_Sinking_Anomalies.SetValue(i,0)
            
        # viscosity.SetValue(i, pviscossity[i]*eta0)
        
        NonDimViscosity.SetValue(i, pviscossity[i])
        
        Composition.InsertTuple4(i,pcompos[i][0],pcompos[i][1],pcompos[i][2],pcompos[i][3])
        
        LLSVp.SetValue(i,pLLSVP[i])
        # stress.InsertTuple6(i, pstress[i][0], pstress[i][1], \
 #                                 pstress[i][2], pstress[i][3], \
 #                                 pstress[i][4], pstress[i][5])
 #
 #        pressure.SetValue(i, ppressure[i]*PressureDim)
 #
 #        composition.SetValue(i, pcomposition[i])
 #
        # radialvelocity.SetValue(i, pvr[i]*velocityDimFactor)
        
        NonDimRadialVelocity.SetValue(i, pvr[i])
        
        radius.SetValue(i, pradius[i]*radiusEarth/1000.0) 
        
        NonDimRadius.SetValue(i, pradius[i]) 
             
        depth.SetValue(i,pdepth[i]*radiusEarth/1000.0)
        
        NonDimDepth.SetValue(i,pdepth[i])
        
        # temperature_Anomaly.SetValue(i, DeltaT*(ptemperature[i]+SurfTemp) - DeltaT*(phoriz_avg[pradius_horiz_avg.index(pradius[i])]+SurfTemp))
        
        NonDimTemperature_Anomaly.SetValue(i, ptemperature[i]-(phoriz_avg[pradius_horiz_avg.index(pradius[i])]))           
                   
        # horizontal_av.SetValue(i, phoriz_avg[pradius_horiz_avg.index(float("{:.5f}".format(pradius[i])))])
        # horizontal_av.SetValue(i, phoriz_avg[pradius_horiz_avg.index(pradius[i])])
        
        # radialheatadv.SetValue(i, ((DeltaT*(ptemperature[i]+SurfTemp)) - (DeltaT*(phoriz_avg[pradius_horiz_avg.index(pradius[i])]+SurfTemp)))*(pvr[i]*velocityDimFactor)/100)
        
        # horizontal_av.SetValue(i,  DeltaT*(phoriz_avg[pradius_horiz_avg.index(pradius[i])]+SurfTemp))
        
        NonDimHorizontal_Av.SetValue(i, (phoriz_avg[pradius_horiz_avg.index(pradius[i])]))
        
        Compositional_Av1.SetValue(i,(pComp_avg_1[pradius_horiz_avg.index(pradius[i])]))
        
        Compositional_Av2.SetValue(i,(pComp_avg_2[pradius_horiz_avg.index(pradius[i])]))
        
        Compositional_Av3.SetValue(i,(pComp_avg_3[pradius_horiz_avg.index(pradius[i])]))
        
        Compositional_Av4.SetValue(i,(pComp_avg_4[pradius_horiz_avg.index(pradius[i])]))
        
        Compositional_Av5.SetValue(i,(pComp_avg_5[pradius_horiz_avg.index(pradius[i])]))
        
        NonDim_Vrms_horiz.SetValue(i,(pVrms_horiz[pradius_horiz_avg.index(pradius[i])]))
        
        NonDim_Vrms_radial.SetValue(i,(pVrms_vertical[pradius_horiz_avg.index(pradius[i])]))
                        
    #end for
    
    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(points)

    # grid.GetPointData().SetVectors(velocity)
        
    # grid.GetPointData().SetActiveScalars('viscosity')
    #     grid.GetPointData().SetScalars(viscosity)
    #
    grid.GetPointData().SetVectors(NonDimVelocity)
    #
    #     grid.GetPointData().SetActiveScalars('temperature')
    #     grid.GetPointData().SetScalars(temperature)

    # grid.GetPointData().SetTensors(stress)
    
    grid.GetPointData().SetActiveScalars('NonDimTemperature')
    grid.GetPointData().SetScalars(NonDimTemperature)
    
    grid.GetPointData().SetActiveScalars('NonDimViscosity')
    grid.GetPointData().SetScalars(NonDimViscosity)
    
    # grid.GetPointData().SetActiveScalars('pressure')
     #    grid.GetPointData().SetScalars(pressure)
     #
    # grid.GetPointData().SetActiveScalars('composition')
     #    grid.GetPointData().SetScalars(composition)
    
    # grid.GetPointData().SetActiveScalars('radialvelocity')
    #     grid.GetPointData().SetScalars(radialvelocity)
    #
    grid.GetPointData().SetActiveScalars('NonDimRadialVelocity')
    grid.GetPointData().SetScalars(NonDimRadialVelocity)
    
    grid.GetPointData().SetActiveScalars('NonDim_Vrms_radial')
    grid.GetPointData().SetScalars(NonDim_Vrms_radial)
    
    grid.GetPointData().SetActiveScalars('NonDim_Vrms_horiz')
    grid.GetPointData().SetScalars(NonDim_Vrms_horiz)
    
    grid.GetPointData().SetActiveScalars('Compositional_Av1')
    grid.GetPointData().SetScalars(Compositional_Av1)
    
    grid.GetPointData().SetActiveScalars('Compositional_Av2')
    grid.GetPointData().SetScalars(Compositional_Av2)
    
    grid.GetPointData().SetActiveScalars('Compositional_Av3')
    grid.GetPointData().SetScalars(Compositional_Av3)
    
    grid.GetPointData().SetActiveScalars('Compositional_Av4')
    grid.GetPointData().SetScalars(Compositional_Av4)
    
    grid.GetPointData().SetActiveScalars('Compositional_Av5')
    grid.GetPointData().SetScalars(Compositional_Av5)
        
    grid.GetPointData().SetActiveScalars('radius')
    grid.GetPointData().SetScalars(radius)
    
    grid.GetPointData().SetActiveScalars('NonDimRadius')
    grid.GetPointData().SetScalars(NonDimRadius)
    
    grid.GetPointData().SetActiveScalars('depth')
    grid.GetPointData().SetScalars(depth)
    
    grid.GetPointData().SetActiveScalars('NonDimDepth')
    grid.GetPointData().SetScalars(NonDimDepth)
    
    # grid.GetPointData().SetActiveScalars('temperature_Anomaly')
  #   grid.GetPointData().SetScalars(temperature_Anomaly)
    
    grid.GetPointData().SetActiveScalars('NonDimTemperature_Anomaly')
    grid.GetPointData().SetScalars(NonDimTemperature_Anomaly)
    
    # grid.GetPointData().SetActiveScalars('radialheatadv')
    #     grid.GetPointData().SetScalars(radialheatadv)
    #
    # grid.GetPointData().SetActiveScalars('horizontal_av')
     #    grid.GetPointData().SetScalars(horizontal_av)
    
    grid.GetPointData().SetActiveScalars('NonDimHorizontal_Av')
    grid.GetPointData().SetScalars(NonDimHorizontal_Av)
    
    # grid.GetPointData().SetActiveScalars('TempAnomalies_5Cases')
    #     grid.GetPointData().SetScalars(TempAnomalies_5Cases)
    
    grid.GetPointData().SetActiveScalars('Hot_Rising_Anomalies')
    grid.GetPointData().SetScalars(Hot_Rising_Anomalies)
    
    grid.GetPointData().SetActiveScalars('Hot_Sinking_Anomalies')
    grid.GetPointData().SetScalars(Hot_Sinking_Anomalies)
    
    grid.GetPointData().SetActiveScalars('Cold_Rising_Anomalies')
    grid.GetPointData().SetScalars(Cold_Rising_Anomalies)
    
    grid.GetPointData().SetActiveScalars('Cold_Sinking_Anomalies')
    grid.GetPointData().SetScalars(Cold_Sinking_Anomalies)
    # vtk.vtkTimeStamp()
    grid.GetPointData().SetActiveScalars('LLSVP')
    grid.GetPointData().SetScalars(LLSVp)
    
    grid.GetPointData().SetActiveScalars('OtherComp')
    grid.GetPointData().SetScalars(Composition)
    
    #connectivity
    for i in range (len(connectivity)):
        hex = vtk.vtkHexahedron()
        for j in range(8):
            hex.GetPointIds().SetId(j, connectivity[i][j])
        #end for
        grid.InsertNextCell(hex.GetCellType(), hex.GetPointIds())
    #end for

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetDataModeToAscii()
    # writer.
    writer.SetInputData(grid) #replaced SetInput with SetInputdata
    
    print(('Writing %s..') % (capFile[34:]+'.vtu'))
    writer.SetFileName(capFile[0:30]+'/Codes_and_Vtks/Vtks_Opt/'+capFile[34:]+'.vtu')
    writer.Write()
    
    
    # hdf5 -->
    # grp = self.__File.create_group("piece%d" % self.__CurrentPiece)
    #     grp.attrs['bounds'] = inp.GetBounds()
    #
    #     grp.create_dataset("cells", data=inp.Cells)
    #     grp.create_dataset("cell_types", data=inp.CellTypes)
    #     grp.create_dataset("cell_locations", data=inp.CellLocations)
    #
    #     pdata = grp.create_group("point_data")
    #     for name in inp.PointData.keys():
    #         pdata.create_dataset(name, data=inp.PointData[name])
    #     s = vtk.vtkRTAnalyticSource()
    #
    #     c = vtk.vtkClipDataSet()
    #     c.SetInputConnection(s.GetOutputPort())
    #     c.SetValue(157)
    #
    #     w = HDF5Writer()
    #     w.SetInputConnection(c.GetOutputPort())
    #     w.SetFileName("test.h5")
    #     w.SetNumberOfPieces(5)
    #
    #     w.Update()
    
    #end

def main():
    
    if(len(sys.argv) != 3):
        print(msg)
        return
    else:
        paramsDict = Core_Citcom.parse_citcoms_cfg_file(sys.argv[2])
        processData(sys.argv[1], paramsDict)
        #end if
    #end if
#end func

if __name__ == "__main__":
    main()

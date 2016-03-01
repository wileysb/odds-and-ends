#-------------------------------------------------------------------------------
# Name:        Topex_model_Adjustment
# Purpose:
#
# Author:      Johannes May, Wiley Bogren
#
# Created:     21/05/2015
# Copyright:   (c) jmay 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Name:        Calc_TOPEX.py
# Purpose:     This program will calculate the vertical angle for a given direction and distance
#  on a given DEM. The maximum angle is calculated for each of 8 directions and summed. A greater angle
#  equates to less exposure.
#
# Author:      olsenk
#
# Created:     17/02/2014
# Copyright:   (c) olsenk 2014
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import arcpy
import os
import csv
import numpy as np
from arcpy import env
from arcpy.sa import *
from math import atan, degrees
arcpy.CheckOutExtension("Spatial")

def buildAzimuthKernel(radiusMeters, referenceDesc): # , ListOut):
    kernel = {}
    for azClassVal in range(0,360,10):
        kernel[azClassVal] = []

    # build kernel from radius and input raster cell size
    window = int(radiusMeters / referenceDesc.meanCellHeight)
    for row in range(-window, window + 1):
        for col in range(-window, window + 1):

            Xdist = abs(col * referenceDesc.meanCellHeight)
            Ydist = abs(row * referenceDesc.meanCellHeight)
            totalDist = pow(Xdist**2 + Ydist**2, 0.5)
            if (totalDist > 0) and (totalDist <= radiusMeters): # and (totalDist % 50 == 0):
                angleDeg = getAzimuth(col,row)
                for azClassVal in range(0,360,10):
                    if (angleDeg >= azClassVal) & (angleDeg < azClassVal + 10):
                        kernel[azClassVal].append([row,col,totalDist])# int(totalDist + 0.5)])

    return(kernel)


def getAzimuth(col,row):
    '''Given a column and row offset (from a center pixel), return the azimuth in degrees.
    0 is north, 90 is east'''
    if row == 0:
        row = 0.000000000001
    angleDeg = int(degrees(atan((-1*col)/float(row))))%360
    if (row > 0) :
        # atan returns 270-90 (north-facing).  if row is negative, the azimuth is actually between 90 and 270 (south-facing)
        angleDeg = (angleDeg + 180)%360

    return angleDeg


def loadPOIlist(poi_fn):
    '''Read in point-type shapefile, get UTM x and y coords, and convert these
    to pixel offsets.'''

    cursor = arcpy.SearchCursor(poi_fn)

    POI_list = []
    for row in cursor:
        x =    row.getValue("Shape").getPart().X
        y =    row.getValue("Shape").getPart().Y
        oid =  row.getValue("FID_")
        POI_list.append([oid, x, y])

    return POI_list


def convertUTMtoPixel(px,py,referenceDesc):
    '''Image coordinates, given image metadata and projected coordinates.

    Assumes coordinates are in the same projection as the image.'''
    originX = referenceDesc.extent.XMin # gt[0]
    dx = referenceDesc.meanCellWidth # gt[1]

    originY  = referenceDesc.extent.YMax # gt[3]
    dy = referenceDesc.meanCellHeight # gt[5]

    x = int((px - originX) / dx) #x pixel
    y = int((py - originY) / (-1*dy)) #y pixel

    return x,y


def Load_window(array_fn, referenceDesc,px,py, radiusMeters = 1000):
    '''Given the array filename and the coordinates of the POI,
    returns a window array around the POI, as well as
    the coordinates of the POI within that window.'''
    windowRadius = int(radiusMeters / referenceDesc.meanCellHeight)

    # find origin
    llx = px - windowRadius
    lly = py - windowRadius

    POIx = windowRadius
    POIy = windowRadius

    windowWidth  = 2*windowRadius + 1
    windowHeight = 2*windowRadius + 1

    # Check for chunks of the window hanging over the image's origin axes
    if llx < 0:
        POIx = POIx + llx
        windowWidth = windowWidth + llx
        llx = 0
    if lly < 0:
        POIy = POIy + lly
        windowHeight = windowHeight + lly
        lly = 0

    # check for chunks of the window hanging over the far edges
    if ( POIx + windowRadius) >referenceDesc.extent.width:
        hangover_x = POIx + windowRadius - referenceDesc.extent.width
        windowWidth = windowWidth - hangover_x
    if ( POIy + windowRadius) >referenceDesc.extent.height:
        hangover_y = POIy + windowRadius - referenceDesc.extent.height
        windowHeight = windowHeight - hangover_y

    array_part = arcpy.RasterToNumPyArray(array_fn, arcpy.Point(llx, lly), windowWidth, windowHeight)

    return array_part, POIx, POIy



def main():
    demImage_fn = "B:\\30-I\\34\\347020 WISLINE\\DEM\\DTED10_Norge_32N_2.tif"
    POI = "B:\\30-I\\34\\347020 WISLINE\\POI\\POI_Z32.shp"

    hdr = ['FID']
    for azClass in range(0,360,10):
        hdr.append('{0} - {1}'.format(azClass,azClass+10))

    out_dir = "B:\\30-I\\34\\347020 WISLINE\\CSVneg"
    out_name = "TopexList_NO_160226.csv"
    out_f = open(os.path.join(out_dir, out_name),'wb')
    out_csv = csv.writer(out_f)
    out_csv.writerow(hdr)

    # read DEM into memory
    referenceDesc = arcpy.Describe(demImage_fn)

    # read POI into memory
    POI_list = loadPOIlist(POI) # list of POI tuples [OID, utm-E, utm-N]

    ###lowerLeftCorner = arcpy.Point(referenceDesc.Extent.XMin + (referenceDesc.meanCellWidth / 2), referenceDesc.Extent.YMin + (referenceDesc.meanCellHeight / 2))
    ###myArray = arcpy.RasterToNumPyArray(demImage, lowerLeftCorner, referenceDesc.width, referenceDesc.height, nodataVal)
    # demImage = arcpy.RasterToNumPyArray(demImage_fn)
    # myArray = arcpy.RasterToNumPyArray(demImage)#, POI, OBJECTID)
    # demImage = np.ma.masked_values(np.rint(myArray), referenceDesc.noDataValue)

    # set workspaces
    arcpy.env.workspace = "c:\\temp"
    arcpy.env.scratchWorkspace = "c:\\temp"

    # pad DEM with nodataVal to avoid boundary issues
##    pad = np.ma.zeros((500,referenceDesc.width), dtype=np.int32)
##    demImage = np.ma.vstack((pad, demImage, pad))
##    referenceDesc.height += 1000
##    pad = np.ma.zeros((referenceDesc.height, 500), dtype=np.int32)
##    demImage = np.ma.hstack((pad, demImage, pad))
##    referenceDesc.width += 1000

    searchKernel = buildAzimuthKernel(1000, referenceDesc)
    x = y = 100 # center coordinates of a 201x201 array
    for POI in POI_list:
        oid = POI[0]
        UTMx = POI[1]
        UTMy = POI[2]

        # Load a window of the DEM with width&height = 1+ 2* radiusPixels = 1+ 2* 100
        # Any portions of the window hanging off the edge of the DEM will be returned as NoData (-3.40282347e+38)
        demImage= arcpy.RasterToNumPyArray(demImage_fn, arcpy.Point(UTMx-1000,UTMy-1000),201,201)
        demImage =  np.where((demImage < -999), np.nan, demImage)

        exposure = {}
        for azClass in range(0,360,10):
            exposure[azClass] = -90
            for pixelOffset in searchKernel[azClass]:
                if ((y + pixelOffset[0]) >= 0 and (y + pixelOffset[0]) < referenceDesc.height and (x + pixelOffset[1]) >= 0 and (x + pixelOffset[1]) < referenceDesc.width):
                    if (int(demImage[y + pixelOffset[0], x + pixelOffset[1]]) != referenceDesc.noDataValue): # Original version tested integer values against noDataValue.  Doing the integer conversion here makes the whole thing a little slower, but therefore we keep our 0.1m height precision
                        exposeArctan = np.degrees(np.arctan((demImage[y + pixelOffset[0],x + pixelOffset[1]] - demImage[y,x]) / (pixelOffset[2]))) # (dZ / dist)
                        exposure[azClass] = max(exposeArctan,exposure[azClass])
        row = [oid]
        for azClass in range(0,360,10):
            row.append(exposure[azClass])
        out_csv.writerow(row)
    out_f.close() # close the output csv




if __name__ == '__main__':
    main()


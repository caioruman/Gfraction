import numpy as np
from osgeo import osr, ogr
import sys

# for reading and writing rpn files
from rpn.rpn import RPN
from rpn import data_types


from netCDF4 import Dataset

"""
Authors: Caio J Ruman

The main function is based on an script written by Sasha.

What the script does:
  Given a list of shapefiles and a domain grid, it reads each feature of each shapefile with the objective of calculating the glacier area and glacier number of each grid cell.

What it needs to run:
  - A rpn file with the corners of the domain (Use the script lolacorners.sh to generate that)
  - A rpn file with the centres of the domain, just so I have the right tic tacs to save the rpn file.
  - A shapefile with the glacier data.

Functions:
 - get_grid_area: Returns the grid area of each grid cell on the domain.
 - get_data_from_shape: Do all the heavy lifting.

List of know problems:
- The script will break if the domain crosses the international date line. It won't handle well a polygon that starts in longitude 179 and ends in the longitude -179. Or lon 1 and Lon 359.
- It is not confirmed, but I suppose that the script will break over the North Pole too.
  -- As there isn't any glacier on those parts, an 'if' will fix the problem.
- The script is inneficient, computationaly speaking. As it loops over the entire grid.
  -- Using a 'seek' function over the latitude/longitude points of the feature and looping throught the nearest grids cell should speed the program considerably.
  -- Done! The runtime went from 10 hours to 10 minutes!!
"""


"""
 feature name list:
   RGIID
   GLIMSID
   RGIFLAG
   BGNDATE
   ENDDATE
   CenLon
   CenLat
   O1Region
   O2Region
   Area
   Zmin
   Zmax
   Zmed
   Slope
   Aspect
   Lmax
   GlacType
   Name
"""

## Input parameters

#input_rpn_file="/glacier/caioruman/Config/DomainSetUp/misc_fields/DX_WC_0.11deg.rpn"
input_rpn_file="/HOME/caioruman/Scripts/LatLonCC/NA022_v2.out.corners"
input_rpn_centres="/HOME/caioruman/Scripts/LatLonCC/NA022_v2.out.centres"
#input_rpn_file="/HOME/caioruman/Scripts/LatLonCC/arcticDomain.out.corners"
#input_rpn_centres="/HOME/caioruman/Scripts/LatLonCC/arcticDomain.out.centres"
#input_shp_file=["RGI_ArcticCanadaNorth_Fixed.shp"]
input_shp_file=["RGI_Greenland_Fixed.shp", "RGI_Alaska_fixed.shp", "RGI_WesternCanadaUS_Fixed.shp", # "RGI_Greenland_Fixed.shp"]
                "RGI_ArcticCanadaNorth_Fixed.shp", "RGI_ArcticCanadaSouth_Fixed.shp"]#, "RGI_Iceland_Fixed.shp",
#                "RGI_CaucasusMiddleEast.shp", "RGI_CentralAsia.shp", "RGI_CentralEurope.shp", "RGI_NorthAsia.shp",
#                "RGI_RussianArctic.shp", "RGI_Scandinavia.shp", "RGI_SouthAsiaEast.shp", "RGI_SouthAsiaWest.shp",
#input_shp_file=["RGI_Svalbard_Fixed.shp"]
#input_shp_file=["/glacier/caioruman/Data/SHP/02_rgi50_WesternCanadaUS.shp", "03_rgi50_ArcticCanadaNorth.shp", "04_rgi50_ArcticCanadaSouth.shp"]#, "01_rgi50_Alaska.shp"]
path_shp = "/glacier/caioruman/Data/SHP/Fixed/"

#output_rpn_file="/glacier/caioruman/Data/GFraction/WC011.rpn"
output_rpn_file="/glacier/caioruman/Data/GFraction/NA022_v2.rpn"

# =============================================

def get_grid_area(lons2d, lats2d):
    """
    Returns the grid area of each grid cell on the domain.
    """
    lons2d_copy = lons2d.copy()
    #lons2d_copy[lons2d > 180] -= 360
    aux = 0
    nx, ny = lons2d.shape

    Grid_area = np.zeros((lons2d.shape[0]-1,lons2d.shape[1]-1), dtype=np.float64)

    for ix in range(nx-1):
        for jy in range(ny-1):
            # Create the geometry for the grid point (i, j), (i+1, j), (i+1, j+1), (i, j+1), (i, j)
            wkt = "POLYGON(({} {}, {} {}, {} {}, {} {}, {} {}))".format(lons2d_copy[ix, jy], lats2d[ix, jy],
                                                                        lons2d_copy[ix+1, jy], lats2d[ix+1, jy],
                                                                        lons2d_copy[ix+1, jy+1], lats2d[ix+1, jy+1],
                                                                        lons2d_copy[ix, jy+1], lats2d[ix, jy+1],
                                                                        lons2d_copy[ix, jy], lats2d[ix, jy])
            #print wkt
            p = ogr.CreateGeometryFromWkt(wkt)

            Grid_area[ix, jy] = p.Area()

    return Grid_area


def get_data_from_shape(lons2d, lats2d, GL_number, GL_area, ZCount, Zmax, Zmin, Zmed, shp_path=""):
    """
    Get the number of Glaciers and area inside a Grid (Geometry) in a shapefile
    """

    ogr.UseExceptions()

    driver = ogr.GetDriverByName("ESRI Shapefile")
    datastore = driver.Open(shp_path, 0)
    layer = datastore.GetLayer(0)

    latlong = osr.SpatialReference()
    latlong.ImportFromProj4("+proj=latlong")

    ct = osr.CoordinateTransformation(latlong, layer.GetSpatialRef())

    print(layer.GetFeatureCount())

    ##read features from the shape file
    feature = layer.GetNextFeature()
    i = 0

    lons2d_copy = lons2d.copy()
    lons2d_copy[lons2d > 180] -= 360
    aux = 0
    nx, ny = lons2d.shape
    # iterate all features in the layer
    while feature:

        # get the geometry of the feature
        geom = feature.GetGeometryRef()
        if i % 250 == 0:
            print "Feature: ", i+1
        #print geom
# put an if here to make the grid skip the north pole


        # Get the center lat and center lons
        clat = feature.GetField("CenLat")
        clon = feature.GetField("CenLon")
        v_zmin = feature.GetField("Zmin")
        v_zmed = feature.GetField("Zmed")
        v_zmax = feature.GetField("Zmax")

        # Getting the point in grid
        a = abs(lats2d-clat)+abs(lons2d_copy-clon)
        ni,nj = np.unravel_index(a.argmin(),a.shape)

        #Checking to see if the point is in the domain. Didn't work. TO DO
        #if (ni < nx) and (nj < ny):
        # Setting up the values of Zmed, Zmin and Zmax
        #    if (Zmin[ni, nj] == 0):
        #        Zmin[ni, nj] = v_zmin
    #        elif v_zmin != 0:
    #            Zmin[ni, nj] = min(v_zmin, Zmin[ni, nj])

    #        Zmax[ni, nj] = max(v_zmax, Zmax[ni, nj])

    #        if (v_zmed != 0):
    #            ZCount[ni, nj] += 1
    #            Zmed[ni, nj] += v_zmed

        """
            # Maybe try to calculate the number of bins based on the slope of each glacier?

        """

        # Now it only loops throught the grids around the first point of the polygon
        #var = geom.ExportToWkt().split(',')[1].split(' ')
        #print geom.ExportToWkt()
        #print var
        #var_lat = float(var[1])
        #var_lon = float(var[0])

        #a = abs(lats2d-var_lat)+abs(lons2d_copy-var_lon)
        #ni,nj = np.unravel_index(a.argmin(),a.shape)

        # Setting the box that I'll loop to put the glacier data

        ni_i = ni - 10
        ni_f = ni + 10
        nj_i = nj - 10
        nj_f = nj + 10

        if (ni_i < 0): ni_i = 0
        if (nj_i < 0): nj_i = 0

        if (ni_f > nx-1): ni_f = nx - 1
        if (nj_f > ny-1): nj_f = ny - 1

        # Looping throught the box

        for ix in range(ni_i, ni_f):

            for jy in range(nj_i, nj_f):
                # I don't want to mess with things in the north pole. But as there isn't glaciers there, this if may be pointless
                if (lats2d[ix, jy] < 88 and lats2d[ix+1, jy] < 88 and lats2d[ix, jy+1] < 88 and lats2d[ix+1, jy+1] < 88):

                    # Create the geometry for the grid point (i, j), (i+1, j), (i+1, j+1), (i, j+1), (i, j)
                    wkt = "POLYGON(({} {}, {} {}, {} {}, {} {}, {} {}))".format(lons2d_copy[ix, jy], lats2d[ix, jy],
                                                                                lons2d_copy[ix+1, jy], lats2d[ix+1, jy],
                                                                                lons2d_copy[ix+1, jy+1], lats2d[ix+1, jy+1],
                                                                                lons2d_copy[ix, jy+1], lats2d[ix, jy+1],
                                                                                lons2d_copy[ix, jy], lats2d[ix, jy])
                    #print wkt
                    try:
                        p = ogr.CreateGeometryFromWkt(wkt)
                    except:
                        #Creating an empty geometry.
                        wkt = "GEOMETRYCOLLECTION EMPTY"
                        p = ogr.CreateGeometryFromWkt(wkt)
                        print "Invalid Grid Box. Check if it is near the international date line"
                        print "#################"

                    # project coordinates to the coordinate system of the shape file
                    #p.Transform(ct)


                    try:
                        intersection = geom.Intersection(p)

                    except:
                        # This will make the polygon valid but will remove the empty areas inside a glacier.
                        # For 1 sample, there is a 8% increase in the glacier area after the "fix"
                        print "Before the fix"
                        print geom
                        geom = geom.Buffer(0)
                        print "#################"
                        print "After the fix"
                        print geom

                        intersection = geom.Intersection(p)

                    if not intersection.IsEmpty():

                        GL_number[ix, jy] += 1
                        GL_area[ix, jy] += intersection.Area()

                        #print intersection
                        #print intersection.Area()
                        #print p.Area()
                        #print "Feature #", i
                        #print "ix, jy:", ix, jy
                        #print "#################################"

                        #print("{}/{}".format(ix + 1, nx))
                        #print("{}/{}".format(jy + 1, ny))
                        #print wkt
                        #print intersection
                        #print geom.Contains(p)

        feature = layer.GetNextFeature()
        i += 1

        #if (i == 50):
        #    break


    datastore.Destroy()
    return GL_number, GL_area, ZCount, Zmax, Zmin, Zmed

def get_lon_lat_fields_from_rpn(rpn_path=""):
    r = RPN(rpn_path)

    var_names = r.get_list_of_varnames()
    sel_name = None
    for vname in var_names:
        if vname not in [">>", "^^", "HY"]:
            sel_name = vname
            break

    r.get_first_record_for_name(sel_name)
    lons, lats = r.get_longitudes_and_latitudes_for_the_last_read_rec()
    r.close()
    return lons, lats

def save_data_to_rpn(files, in_file="", out_file=""):

    rin = RPN(in_file)
    rout = RPN(out_file, mode="w")


    # Read coordinates and reshape(needed for writing)
    x = rin.get_first_record_for_name(">>")
    x.shape = (-1, 1)
    print(x.shape)


    y = rin.get_first_record_for_name("^^")
    y.shape = (1, -1)

    # get parameters of the last read record
    coord_info = rin.get_current_info()

    print(coord_info)

    # write coordinates
    coord_info.update({"name": ">>", "label": "NGP", "typ_var": "X", "nbits": -coord_info["nbits"]})
    rout.write_2d_field_clean(x, properties=coord_info)

    coord_info.update({"name": "^^"})
    rout.write_2d_field_clean(y, properties=coord_info)

    # write the data
    #for item in files:
    #rout.write_2d_field_clean(files, properties=dict(name="DX", label="Area", ig=coord_info["ip"] + [0,]))

    for item in files:
        rout.write_2d_field_clean(item[0], properties=dict(name=item[1], label=item[2], ig=coord_info["ip"] + [0,]))

    rin.close()
    rout.close()


# The main function
def main():
    lons2d, lats2d = get_lon_lat_fields_from_rpn(rpn_path=input_rpn_file)

    GL_number = np.zeros((lons2d.shape[0]-1,lons2d.shape[1]-1), dtype=np.int)
    ZCount = np.zeros((lons2d.shape[0]-1,lons2d.shape[1]-1), dtype=np.int)
    Zmax = np.zeros((lons2d.shape[0]-1,lons2d.shape[1]-1), dtype=np.float64)
    Zmin = np.zeros((lons2d.shape[0]-1,lons2d.shape[1]-1), dtype=np.float64)
    Zmed = np.zeros((lons2d.shape[0]-1,lons2d.shape[1]-1), dtype=np.float64)
    GL_area = np.zeros((lons2d.shape[0]-1,lons2d.shape[1]-1), dtype=np.float64)

    Grid_area = get_grid_area(lons2d, lats2d)

    #save_data_to_rpn(Grid_area, in_file=input_rpn_centres, out_file=output_rpn_file)

    for item in input_shp_file:
        print "{}{}".format(path_shp, item)
        GL_number, GL_area, ZCount, Zmax, Zmin, Zmed = get_data_from_shape(lons2d, lats2d, GL_number, GL_area, ZCount, Zmax, Zmin, Zmed, shp_path="{}{}".format(path_shp, item))

    ZCount[ZCount == 0] = 1
    Zmed = Zmed / ZCount

    items = [(GL_number, "GLNM", "GL_number"), (GL_area, "GAREA", "GL_area"), (Grid_area, "DXDY", "Grid_area"),
             (GL_area/Grid_area, "GLF", "GL_fraction"), (Zmed, "ZMED", "Z_Average"), (Zmax, "ZMAX", "Z_Max"),
             (Zmin, "ZMIN", "Z_Min")]

    save_data_to_rpn(items, in_file=input_rpn_centres, out_file=output_rpn_file)


# Entry point of the script
if __name__ == "__main__":
    main()

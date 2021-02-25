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
input_rpn_file="/HOME/caioruman/Scripts/LatLonCC/BC_Alberta.out.corners"
#input_rpn_centres="/HOME/caioruman/Scripts/LatLonCC/NA022_v2.out.centres"
#input_rpn_file="/HOME/caioruman/Scripts/LatLonCC/fileStep1.NA044.out.corners"
#input_rpn_file="/HOME/caioruman/Scripts/LatLonCC/WC011.out.corners"
#input_rpn_centres="/HOME/caioruman/Scripts/LatLonCC/fileStep1.NA044.out.centres"
#input_rpn_centres="/HOME/caioruman/Scripts/LatLonCC/WC011.out.centres"
input_rpn_centres="/HOME/caioruman/Scripts/LatLonCC/BC_Alberta.out.centres"
#input_rpn_file="/HOME/caioruman/Scripts/LatLonCC/arcticDomain.out.corners"
#input_rpn_centres="/HOME/caioruman/Scripts/LatLonCC/arcticDomain.out.centres"
#input_shp_file=["RGI_ArcticCanadaNorth_Fixed.shp", "RGI_ArcticCanadaSouth_Fixed.shp"]
input_shp_file=["RGI_Alaska_fixed.shp", "RGI_WesternCanadaUS_Fixed.shp"]# "RGI_Greenland_Fixed.shp"
#                "RGI_ArcticCanadaNorth_Fixed.shp", "RGI_ArcticCanadaSouth_Fixed.shp"]
#               "RGI_Iceland_Fixed.shp", "RGI_Svalbard_Fixed.shp","RGI_RussianArctic.shp", "RGI_NorthAsia.shp",
#                "RGI_CaucasusMiddleEast.shp", "RGI_CentralAsia.shp", "RGI_CentralEurope.shp", "RGI_NorthAsia.shp",
#                "RGI_RussianArctic.shp", "RGI_Scandinavia.shp", "RGI_SouthAsiaEast.shp", "RGI_SouthAsiaWest.shp",
#input_shp_file=["RGI_Svalbard_Fixed.shp"]
#input_shp_file=["/glacier/caioruman/Data/SHP/02_rgi50_WesternCanadaUS.shp", "03_rgi50_ArcticCanadaNorth.shp", "04_rgi50_ArcticCanadaSouth.shp"]#, "01_rgi50_Alaska.shp"]

#c_values = [0.31, 0.31, 0.37, 0.31, 0.31]
c_values = [0.2205, 0.2205, 0.2205, 0.2205, 0.2205, 0.2205, 0.2205, 0.2205]

path_shp = "/glacier/caioruman/Data/SHP/Fixed/"

#output_rpn_file="/glacier/caioruman/Data/GFraction/WC011.rpn"
#output_rpn_file="/glacier/caioruman/Data/GFraction/NA022_v3_1.375_newc.rpn"
output_rpn_file="/glacier/caioruman/Data/GFraction/BC_Alberta.out.RGI.rpn"

fileout = '/glacier/caioruman/Data/GFraction/BC_Alberta.out.RGI.dat'
saida = open(fileout, 'w')
saida.write("Glacier Volume by region\n")
saida.close()

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


def get_data_from_shape(lons2d, lats2d, GL_number, GL_area, GL_areaReal, GL_volume, Grid_areaReal, Grid_area, GL_Height, c_value, shp_path):
    """
    Get the number of Glaciers and area inside a Grid (Geometry) in a shapefile
    """

    ogr.UseExceptions()

    volumeTotal = 0
    areaTotal = 0

    saida = open(fileout, 'a')
    saida.write(shp_path + "\n")

    driver = ogr.GetDriverByName("ESRI Shapefile")
    datastore = driver.Open(shp_path, 0)
    layer = datastore.GetLayer(0)

    latlong = osr.SpatialReference()
    latlong.ImportFromProj4("+proj=latlong")

    ct = osr.CoordinateTransformation(latlong, layer.GetSpatialRef())

    #gamma = 1.357
    gamma = 1.375
    gamma_sheets = 1.25
    c = 0.2055
    c_sheets = 1.7026

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

        # Get the center lat and center lons
        clat = feature.GetField("CenLat")
        clon = feature.GetField("CenLon")
        # Area comes in squared Km, must transform to squared m
        area = feature.GetField("Area")*1000000
        v_zmin = feature.GetField("Zmin")
        v_zmed = feature.GetField("Zmed")
        v_zmax = feature.GetField("Zmax")
        # Area in unknow units
        areaGeom = geom.Area()

        # If a glacier is greater than the grid area, its an ice cap
        # See the paper and change this value based on the area
        if (area > 1000000000):
            volume = c_sheets * area**gamma_sheets
        else:
            volume = c_value * area**gamma

        volumeTotal += volume
        areaTotal += area

        ratio = area/areaGeom

#        print "Area feature:", area
#        print "Area Geom", areaGeom
#        print "ratio", area/areaGeom
#        print "Volume feature", c * (area)**gamma
#        print "Volume python", c * (areaGeom*ratio)**gamma

        # Getting the point in grid
        a = abs(lats2d-clat)+abs(lons2d_copy-clon)
        ni,nj = np.unravel_index(a.argmin(),a.shape)

        """
            # Maybe try to calculate the number of bins based on the slope of each glacier?

        """
        print clat, clon
        print a
        sys.exit()
        # Setting the box that I'll loop to put the glacier data

        # *4 for 0.11, *2 for 0.22, *1 for 0.44
        ni_i = ni - 10*2
        ni_f = ni + 10*2
        nj_i = nj - 10*2
        nj_f = nj + 10*2

#        ni_i = ni - 1
#        ni_f = ni + 1
#        nj_i = nj - 1
#        nj_f = nj + 1

        if (ni_i < 0): ni_i = 0
        if (nj_i < 0): nj_i = 0

        if (ni_f > nx-1): ni_f = nx - 1
        if (nj_f > ny-1): nj_f = ny - 1

        # Looping throught the box

        area_test = 0

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
                        GL_areaReal[ix, jy] += intersection.Area()*ratio
                        GL_Height[ix, jy] += v_zmed

                        area_test += intersection.Area()*ratio

                        # if the glacier area is almost the same site as the grid area
                        #if (intersection.Area() >= Grid_area[ix, jy]*0.95):
                            # Ice cap or Ice Sheet
                        #    GL_volume[ix, jy] += c * (intersection.Area()*ratio)**gamma_sheets
                        #else:
                            # Single Glaciers
                        GL_volume[ix, jy] += intersection.Area()*volume/areaGeom

                        Grid_areaReal[ix, jy] = Grid_area[ix, jy]*ratio

                        #print "Area real", intersection.Area()*ratio
                        #print "Area python", intersection.Area()
                        #print "Volume real", GL_volume[ix, jy]
                        #print "Area Grid real", Grid_area[ix, jy]*ratio
                        #print "Area grid python", Grid_area[ix, jy]
                        #print "i", ix,  "j", jy



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

#        if not (-1 < (area - area_test) < 1):
#            print area - area_test
#            print "Area from feature: ", area
#            print "Area put in the grids:", area_test

        feature = layer.GetNextFeature()
        volume_aux = 0
        #sys.exit(0)
        i += 1

        #if (i == 50):
        #    break


    datastore.Destroy()

    print volumeTotal
    saida.write("Volume (m3): {0}\n".format(volumeTotal))
    saida.write("Area (m2): {0}\n".format(areaTotal))
    saida.close()
    return GL_number, GL_area, GL_areaReal, GL_volume, Grid_areaReal, GL_Height

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
    coord_info.update({"varname": ">>", "label": "NGP", "typ_var": "X", "nbits": -coord_info["nbits"]})
    rout.write_2d_field_clean(x, properties=coord_info)

    coord_info.update({"varname": "^^"})
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
    GL_areaReal = np.zeros((lons2d.shape[0]-1,lons2d.shape[1]-1), dtype=np.float64)
    GL_volume = np.zeros((lons2d.shape[0]-1,lons2d.shape[1]-1), dtype=np.float64)
    GL_area = np.zeros((lons2d.shape[0]-1,lons2d.shape[1]-1), dtype=np.float64)
    Grid_areaReal = np.zeros((lons2d.shape[0]-1,lons2d.shape[1]-1), dtype=np.float64)
    GL_Height = np.zeros((lons2d.shape[0]-1,lons2d.shape[1]-1), dtype=np.float64)

    Grid_area = get_grid_area(lons2d, lats2d)

    #save_data_to_rpn(Grid_area, in_file=input_rpn_centres, out_file=output_rpn_file)

    for item, c_value in zip(input_shp_file, c_values):
        print "{}{}".format(path_shp, item)
        GL_number, GL_area, GL_areaReal, GL_volume, Grid_areaReal, GL_Height = get_data_from_shape(lons2d, lats2d, GL_number, GL_area, GL_areaReal, GL_volume, Grid_areaReal, Grid_area, GL_Height, c_value, shp_path="{}{}".format(path_shp, item))

    #ZCount[ZCount == 0] = 1
    #Zmed = Zmed / ZCount

    GL_number_height = GL_number.copy()
    GL_number_height[GL_number_height == 0] = 1

    items = [(GL_number, "GLNM", "GL_number"), (GL_area, "GAREA", "GL_area"), (Grid_area, "DXDY", "Grid_area"),
             (GL_area/Grid_area, "GLF", "GL_fraction"), (GL_areaReal, "GLA", "GL_Area"), (GL_volume, "GVOL", "GL_Volume"),
             (GL_areaReal, "DX", "Grid_areaR"), (GL_Height/GL_number_height, "GLH", "GL_MHeigh")]

    save_data_to_rpn(items, in_file=input_rpn_centres, out_file=output_rpn_file)




# Entry point of the script
if __name__ == "__main__":
    main()

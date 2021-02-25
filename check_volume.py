import numpy as np
from osgeo import osr, ogr
import sys

# for reading and writing rpn files
from rpn.rpn import RPN
from rpn import data_types


from netCDF4 import Dataset

"""
Authors: Caio J Ruman


Estimate the Glacier Volume over greenland.
"""

## Input parameters

#input_rpn_file="/glacier/caioruman/Config/DomainSetUp/misc_fields/DX_WC_0.11deg.rpn"
input_rpn_file1="/glacier/caioruman/Data/Geophys/NA022/geophysNA022_CRUv4"


output_rpn_file="/glacier/caioruman/Data/GFraction/NA022.Caiov2.greenland"
input_rpn_file2="/glacier/caioruman/Data/GFraction/NA022.Caiov2.rpn"

# =============================================



def calculate_gvol(lons2d, lats2d, GVOL, input_rpn_file1, input_rpn_file2):
    """
    Check for the glacier volume
    """

    # GLNM
    r = RPN(input_rpn_file1)
    #glnm = r.get_first_record_for_name("GLNM")
    glf = r.get_4d_field("VF")
    dates = list(glf.keys())
    #getting only the MASK values, removing the date dict
    glf = glf[dates[0]]
    values = list(glf.keys())

    # putting the data in a np.array
    glf = np.array([glf[values[x]] for x in range(0,len(values),1)])
    print glf.shape
    glf = glf[1]



#    glnm = r.get_4d_field("GLNM")
#    dates = list(glnm.keys())
#    #getting only the MASK values, removing the date dict
#    glnm = glnm[dates[0]]
#    values = list(glnm.keys())
#    # putting the data in a np.array
#    glnm = np.array([glnm[values[x]] for x in range(0,len(values),1)])[0]

#    print "Glacier Number read"
    dxdy = r.get_4d_field("DXDY")
    dates = list(dxdy.keys())
    #getting only the MASK values, removing the date dict
    dxdy = dxdy[dates[0]]
    values = list(dxdy.keys())
    # putting the data in a np.array
    dxdy = np.array([dxdy[values[x]] for x in range(0,len(values),1)])[0]

    print "Grid Area Read"
    # complete the code here
    # GVOL
    r = RPN(input_rpn_file2)
    #gvol_file = r.get_first_record_for_name("GVOL")
    gvol_file = r.get_4d_field("GVOL")
    dates = list(gvol_file.keys())
    #getting only the MASK values, removing the date dict
    gvol_file = gvol_file[dates[0]]
    values = list(gvol_file.keys())
    # putting the data in a np.array
    gvol_file = np.array([gvol_file[values[x]] for x in range(0,len(values),1)])[0]


#    r = RPN('/glacier/caioruman/Data/Geophys/NA022/mask_greenland')
#    mask = r.get_4d_field("MASK")
#    dates = list(mask.keys())
#    #getting only the MASK values, removing the date dict
#    mask = mask[dates[0]]
#    values = list(mask.keys())
#    # putting the data in a np.array
#    mask = np.array([mask[values[x]] for x in range(0,len(values),1)])[0]

    print "Glacier Volume read"

    gamma = 1.25
    gamma_1 = 1.375
    c = 1.7026*1.5
    c_1 = 0.2205*.65

    aux = 0
    nx, ny = lons2d.shape
    print lons2d.shape
    print gvol_file.shape
    print dxdy.shape
    # check all the points
    print xrange(0,nx)
    print xrange(0, ny)

    for i in xrange(0,nx):
        for j in xrange(0, ny):
            if (glf[i,j] > 0) and (gvol_file[i,j] == 0):
                if (glf[i,j] > 0.96):
                    GVOL[i,j] = c * (dxdy[i,j]*glf[i,j])**gamma
                else:
                    GVOL[i,j] = c_1 * (dxdy[i,j]*glf[i,j])**gamma_1
#            !if (mask[i,j] != 0) and (glnm[i,j] > 0):
#            !    if (glf[i,j] >= 0.99):
#            !        GVOL[i,j] = c * dxdy[i, j]**gamma
#            !    else:
#            !        GVOL[i,j] = c_1 * (dxdy[i, j]*glf[i,j])**gamma_1
#            !    if (gvol_file[i,j] != 0) and (glf[i,j] < 0.9):
#            !        print gvol_file[i,j], GVOL[i,j]
#            !        print "ratio", gvol_file[i,j]/GVOL[i,j]
#            !        if gvol_file[i,j]/GVOL[i,j] > 0.1:
#            !            GVOL[i,j] = gvol_file[i,j]
#            !if (mask[i,j] == 0):
#            !    GVOL[i,j] = gvol_file[i,j]

    return GVOL+gvol_file

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
    lons2d, lats2d = get_lon_lat_fields_from_rpn(rpn_path=input_rpn_file2)

    GVOL = np.zeros((lons2d.shape[0],lons2d.shape[1]), dtype=np.float64)


    GVOL = calculate_gvol(lons2d, lats2d, GVOL, input_rpn_file1, input_rpn_file2)


    items = [(GVOL, "GVOL", "GL_Volume")]

    save_data_to_rpn(items, in_file=input_rpn_file2, out_file=output_rpn_file)


# Entry point of the script
if __name__ == "__main__":
    main()

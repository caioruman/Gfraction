# Gfraction
Given a list of shapefiles and a domain grid, it reads each feature of each shapefile with the objective of calculating the glacier area and glacier number of each grid cell.

What it needs to run:
  - A rpn file with the corners of the domain (Use the script lolacorners.sh to generate that)
  - A rpn file with the centres of the domain, just so I have the right tic tacs to save the rpn file.
  - A shapefile with the glacier data.

Functions:
 - get_grid_area: Returns the grid area of each grid cell on the domain.
 - get_data_from_shape: Do all the heavy lifting.

Code from Jan-2017

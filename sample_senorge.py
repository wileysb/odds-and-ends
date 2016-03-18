import argparse
import ogr
import osr
import csv
import netCDF4
import numpy as np
import datetime as dt

utm33n = osr.SpatialReference()
utm33n.ImportFromEPSG(32633)

# Define loading string for senorge files
metv2 = '/space/wib_data/CLIMATE/METNO/senorgeV2/temperature/' \
        '{0}{1}/seNorge_v2_0_TEMP1d_grid_{0}{1}{2}.nc'  # .format(year,str(month).zfill(2),str(day).zfill(2)

# Load one file, get X and Y coordinates (UTM33)
tam1 = netCDF4.Dataset(metv2.format(1999, '08', '02'), 'r')
xc = tam1.variables['X'][:] - 500  # left cell edge??
yc = tam1.variables['Y'][:] + 500  # Top cell edge??


def define_command_parser():
    parser = argparse.ArgumentParser()
    parser.description = 'This is a utility for reading a csv file of point coordinates, ' \
                         'and writing a file for seNorge values matching the specified coordinates and time range.\n' \
                         'Unfortunately, the seNorge grids are hardcoded into the utility.'
    parser.add_argument('-i', '--inputfile', type=str, help='Input coordinate list', required=True)
    parser.add_argument('-o', '--outputfile', type=str, help='Output file', required=True)
    parser.add_argument('-f', '--firstyear', type=int, help='First year in timerange to extract')
    parser.add_argument('-l', '--lastyear', type=int, help='Last year in timerange to extract')
    parser.add_argument('-x', '--xcol', type=int, help='Column number of X coordinates')
    parser.add_argument('-y', '--ycol', type=int, help='Column number of Y coordinates')
    parser.add_argument('-u', '--utmcol', type=int, help='Column to read for UTM zone, if varying.  Assumes North')
    parser.add_argument('-E', '--EPSG', type=int, help='EPSG of coordinate system, if fixed; default 4326 (Geographic)')
    # parser.add_argument('-e', '--elevation', type=bool, help='Whether to add senorge cell elevation to the output')
    parser.add_argument('-d', '--dset', type=str, help='temp or precip?')
    return parser


def sort_coordinates(input_sites, args):
    for site in input_sites:
        if args.EPSG == 32633:  # All coordinates are already in UTM33N
            input2utm33 = None
        else:
            if args.utmcol:  # Coordinates are in UTM, but the zone varies by column
                in_epsg = int('326'+str(site[args.utmcol-1]))
                if in_epsg == 32633:  # This row's coordinates are already in UTM33N
                    input2utm33 = None
                else:  # This row's coordinates are in a different zone; define a transform
                    input2utm33 = return_transform(in_epsg)
            else:  # Coordinates are not in UTM; define a transform
                input2utm33 = return_transform(args.EPSG)
        prepend_to_site(site, args, input2utm33)
    return input_sites


def return_transform(in_epsg):
    input_coordinate_system = osr.SpatialReference()
    input_coordinate_system.ImportFromEPSG(in_epsg)
    input2utm33 = osr.CoordinateTransformation(input_coordinate_system, utm33n)
    return input2utm33


def prepend_to_site(site, args, input2utm33):
    in_x = float(site[args.xcol-1])
    in_y = float(site[args.ycol-1])
    if input2utm33:
        point = ogr.CreateGeometryFromWkt("Point ({0} {1})".format(in_x, in_y))  # x,y
        point.Transform(input2utm33)
        in_y = round(point.GetY(), 1)
        in_x = round(point.GetX(), 1)
    site.insert(0, in_y)
    site.insert(0, in_x)
    x_off = np.where(site[0] < xc)[0][0]-1
    y_off = np.where(site[1] < yc)[0][-1]+1
    site.insert(0, y_off)
    site.insert(0, x_off)


if __name__ == '__main__':

    # Parse commandline arguments
    argparser = define_command_parser()
    args = argparser.parse_args()

    # Elevation!

    # Temperature/Precipitation?
    varname = 'mean_temperature'  # or whatever it is for precipitation

    # Load input coordinates
    input_file = open(args.inputfile, 'r')
    input_reader = csv.reader(input_file, dialect='excel')

    # Split into column headers, and data
    input_header = input_reader.next()
    input_sites = [site for site in input_reader]

    input_file.close()
    # Translate coordinates to UTM33, and x/y pixel offsets
    # input_header.insert(0, ' ')  # First column
    input_header.insert(0, 'utm33n_y')
    input_header.insert(0, 'utm33n_x')
    input_header.insert(0, 'y_off')
    input_header.insert(0, 'x_off')

    input_sites = sort_coordinates(input_sites, args)

    # Write output, starting with the input columns.
    output_file = open(args.outputfile, 'wb')
    writer = csv.writer(output_file, dialect='excel')
    header_block = map(list, zip(*input_sites))  # Output format is transposed; columns -> rows
    for i in range(4, len(input_header)):
        header_block[i].insert(0, input_header[i])
        writer.writerow(header_block[i])

    # Iterate over dates
    year = args.firstyear
    i_date = dt.date(args.firstyear, 1, 1)
    while year < args.lastyear+1:
        # Set date variables from iterating date
        year = i_date.year
        month = str(i_date.month).zfill(2)
        day = str(i_date.day).zfill(2)

        # Load netcdf
        nc_day = netCDF4.Dataset(metv2.format(year, month, day), 'r')
        out = [i_date.isoformat(), ]

        # Iterate over seed_sites
        for site in range(len(input_sites)):
            site_day_val = nc_day.variables[varname][0, input_sites[site][1], input_sites[site][0]]
            out.append(str(round(float(site_day_val), 1)))

        # Print mean temperature at each site for this date
        writer.writerow(','.join(out))

        # Iterate date
        i_date = i_date+dt.timedelta(days=1)
        year = i_date.year

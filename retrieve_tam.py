import ogr
import osr
import netCDF4
import numpy as np
import datetime as dt

seeds_fn = '/home/sl_wib/seed_climate/seed_sites.csv'
seeds = np.loadtxt(seeds_fn,skiprows=1,delimiter=',')

thredds = 'http://thredds.met.no/thredds/dodsC/arcticdata/met.no/tam24hNOgrd1957on/tam24hNOgrd1957on_{0}_{1}_{2}.nc'#.format(year,str(month).zfill(2),str(day).zfill(2)
# tam1 = netCDF4.Dataset(thredds.format(1999,'08','02'),'r')
# xc = tam1.variables['Xc'][:] - 500 # left cell edge??
# yc = tam1.variables['Yc'][:] + 500 # Top cell edge??

metv2 = '/space/wib_data/CLIMATE/METNO/senorgeV2/temperature/{0}{1}/seNorge_v2_0_TEMP1d_grid_{0}{1}{2}.nc'#.format(year,str(month).zfill(2),str(day).zfill(2)
tam1 = netCDF4.Dataset(metv2.format(1999,'08','02'),'r')
xc = tam1.variables['X'][:] - 500 # left cell edge??
yc = tam1.variables['Y'][:] + 500 # Top cell edge??

### convert decimal degree coordinates to utm33
utm33n = osr.SpatialReference()
utm33n.ImportFromEPSG(32633)

wgs84 = osr.SpatialReference()
wgs84.ImportFromEPSG(4326)

wgs2utm33 = osr.CoordinateTransformation(wgs84,utm33n)

def ll2utm33(lon, lat):
    """utm_e,utm_n = ll2utm33(lon,lat)"""
    point = ogr.CreateGeometryFromWkt("Point ({0} {1})".format(lon,lat))
    point.Transform(wgs2utm33)
    return point.GetX(),point.GetY()

# convert seed coords to utm33
site_utms = seeds.copy()
for i in range(len(seeds)):
    site_utms[i,1:] = ll2utm33(seeds[i,1],seeds[i,2])


# convert utm seed coords to senorge pixel offsets (tam[x,y])
site_offs = site_utms.copy()
for i in range(len(seeds)):
    x_off = np.where(site_utms[i,1]<xc)[0][0]-1
    y_off = np.where(site_utms[i,2]<yc)[0][-1]+1 
    site_offs[i,1:] = x_off,y_off


# Make header ['Date', site1, site2, site3...]
hdr = site_offs[:,0] 
print ','.join(['Date',]+[str(int(x)) for x in hdr])

# Define date range
start_year = 1960
end_year   = 2006

# Iterate over dates
year=start_year
i_date = dt.date(start_year,1,1)
while year<end_year+1:
    # Set date variables from iterating date
    year = i_date.year
    month = str(i_date.month).zfill(2)
    day   = str(i_date.day).zfill(2)

    # Load netcdf
    nc_day = netCDF4.Dataset(metv2.format(year,month,day),'r')
    out = [i_date.isoformat(),]

    # Iterate over seed_sites
    for seed_site in range(len(hdr)):
        out.append(nc_day.variables['tam'][site_offs[seed_site,2],site_offs[seed_site,1]])

    # Print mean temperature at each site for this date
    print ','.join(out)
    # Iterate date
    i_date=i_date+dt.timedelta(days=1)


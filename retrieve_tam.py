import ogr
import osr
import netCDF4
import os
import numpy as np
import datetime as dt

seeds_fn = '/home/sl_wib/seed_climate/seed_sites.csv'
seeds = np.loadtxt(seeds_fn,skiprows=1,delimiter=',')

thredds = 'http://thredds.met.no/thredds/dodsC/arcticdata/met.no/tam24hNOgrd1957on/tam24hNOgrd1957on_{0}_{1}_{2}.nc'#.format(year,str(month).zfill(2),str(day).zfill(2)
metv2 = '/space/wib_data/CLIMATE/METNO/senorgeV2/temperature/{0}{1}/seNorge_v2_0_TEMP1d_grid_{0}{1}{2}'#.format(year,str(month).zfill(2),str(day).zfill(2)
# tam_date = netCDF4.Dataset(thredds.format(year,str(month).zfill(2),str(day).zfill(2)),'r')
# tam = tam_date.variables['tam'][y,x] # or [x,y]?
# tam1 = netCDF4.Dataset(thredds.format(1999,'08','02'),'r')
# xc = tam1.variables['Xc'][:] - 500 # left cell edge??
# yc = tam1.variables['Yc'][:] + 500 # Top cell edge??
tam1 = netCDF4.Dataset(metv2.format(1999,'08','02'),'r')
xc = tam1.variables['X'][:] - 500 # left cell edge??
yc = tam1.variables['Y'][:] + 500 # Top cell edge??

### convert decimal degree coordinates to utm33
utm33n = osr.SpatialReference()
utm33n.ImportFromEPSG(32633)

wgs84 = osr.SpatialReference()
wgs84.ImportFromEPSG(4326)

wgs2utm33 = osr.CoordinateTransformation(wgs84,utm33n)

def ll2utm33(lon,lat):
    '''utm_e,utm_n = ll2utm33(lon,lat)'''
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

# Loop through dates!
hdr = site_offs[:,0] 
print ','.join(['Date',]+[str(x) for x in hdr])
start_year = 1960
end_year   = 2006

#year=start_year
#num_dates=0
#i_date = dt.date(start_year,1,1)
#while year<end_year+1:
#    i_date=i_date+dt.timedelta(days=1)
#    num_dates+=1
#    year=i_date.year

# make output array

year=start_year
i_date = dt.date(start_year,1,1)
while year<end_year+1:
    # hdr += i_date.isoformat()
    year = i_date.year
    month = str(i_date.month).zfill(2)
    day   = str(i_date.day).zfill(2)
    i_date=i_date+dt.timedelta(days=1)
    nc_day = netCDF4.Dataset(thredds.format(year,month,day),'r')
    out = []
    for seed_site in range(len(hdr)):
        # out_array[seed_site,len(hdr)] = nc_day.variables['tam'][site_offs[seed_site,2],site_offs[seed_site,1]]
        out.append(
# for each day in 1960-2006:
#     nc_day = netCDF4.Dataset(thredds.format(year,str(month).zfill(2),str(day).zfill(2)),'r')
#     tam_day = nc_day.variables['tam'][y,x]
#     for each seed site:
#         add the day's mean temp to the site's timeseries
# Write out the site temp series

import ogr
import osr
import sys

buffer_distance = int(sys.argv[1])

#buffers = [20,50,100,200,500,1000]
tile_index_ds = '/space/wib_data/LANDCOVER/tile_overview/tile_bounds.shp'
pts_ds = '/home/sl_wib/ss_ar_buffers/utm_koor.shp'
out_shp = '/home/sl_wib/ss_ar_buffers/buffer_out/buffer{0}m.shp'#.format(buffer_distance)
ss_ar_fmt = '/space/wib_data/LANDCOVER/ss_ar_shp/ss_ar_{}.shp'#.format(tile)
utm33 = osr.SpatialReference()
utm33.ImportFromEPSG(32633)

# Load points
shpdr = ogr.GetDriverByName('ESRI Shapefile')
pts_src = shpdr.Open(pts_ds,0)
pts_lyr = pts_src.GetLayer()

# Load tile index
tile_src = shpdr.Open(tile_index_ds,0)
tile_idx_lyr = tile_src.GetLayer()


# Build some functions
def Mk_outfile(dst_shp):
    tile_src = shpdr.Open(ss_ar_fmt.format(4599), 0)
    tile_lyr = tile_src.GetLayer()

    dst = shpdr.CreateDataSource(dst_shp)
    dst_lyr = dst.CreateLayer('',None,ogr.wkbPolygon)
    src_lyr_defn = tile_lyr.GetLayerDefn()
    num_fields = src_lyr_defn.GetFieldCount()
    field_names = []
    for field_i in range(num_fields):
        src_field_defn = src_lyr_defn.GetFieldDefn(field_i)
        field_name = src_field_defn.GetName()
        if field_name not in ['cell_id','feature_id','ctr_x_cell','ctr_y_cell','ctr_x_poly','ctr_y_poly','area_m2']:
            field_names.append(field_name)
            dst_lyr.CreateField(src_field_defn)
    dst_lyr.CreateField(ogr.FieldDefn('area',ogr.OFTReal))
    dst_lyr.CreateField(ogr.FieldDefn('PlotID',ogr.OFTString))
    dst_lyr.CreateField(ogr.FieldDefn('BuffDist',ogr.OFTInteger))
    dst_lyr_defn = dst_lyr.GetLayerDefn()

    return dst, dst_lyr, dst_lyr_defn, field_names


def Get_intersecting_tiles(buff_circ,tile_idx_lyr):
    # tile = tile_idx_lyr.GetFeature(0)
    intersecting_tiles=[]
    # while tile:
    for tile_i in range(tile_idx_lyr.GetFeatureCount()):
        tile = tile_idx_lyr.GetFeature(tile_i)
        tile_geom = tile.GetGeometryRef()
        if tile_geom.Intersects(buff_circ):
            intersecting_tiles.append(tile.GetField('tile_id'))
        # tile=tile_idx_lyr.GetNextFeature()
    return intersecting_tiles


def Stamp_one_cookie(pt_feature, tile_idx_lyr, dst_lyr, dst_lyr_defn, field_names, buffer_distance):
    buf_ctr = pt_feature.GetGeometryRef()
    buff_circ = buf_ctr.Buffer(buffer_distance)
    intersecting_tiles = Get_intersecting_tiles(buff_circ,tile_idx_lyr)
    if len(intersecting_tiles) is not 0:
        for tile in intersecting_tiles:
            tile_src = shpdr.Open(ss_ar_fmt.format(tile),0)
            tile_lyr = tile_src.GetLayer()
            for ss_ar_fid in range(0,tile_lyr.GetFeatureCount()):
                src_feat = tile_lyr.GetFeature(ss_ar_fid)
                ss_ar_geom = src_feat.GetGeometryRef()
                if ss_ar_geom.Intersects(buff_circ):
                    dest_geom = ss_ar_geom.Intersection(buff_circ)
                    dest_feat = ogr.Feature(dst_lyr_defn)
                    dest_feat.SetGeometry(dest_geom)
                    dest_feat.SetField('area',dest_geom.GetArea())
                    dest_feat.SetField('PlotID',pt_feature.GetField(0))
                    dest_feat.SetField('BuffDist',buffer_distance)
                    for field_name in field_names:
                        field_val = src_feat.GetField(field_name)
                        dest_feat.SetField(field_name, field_val)
                    dst_lyr.CreateFeature(dest_feat)
    return dst_lyr


def Mk_proj(prj,outf):
    '''Write an ESRI style .prj file defining
    spatial reference for a shapefile dataset.
    :param prj: osr.SpatialReference() object set to dataset projection
    :param outf: (str) path to shapefile dataset, no extensions
    :return: None
    '''
    with open("%s.prj" % outf, "w") as proj:
        prj_out = prj.ExportToPrettyWkt()
        proj.write(prj_out)


#for buf_i in range(len(buffers)):
#    buffer_distance = buffers[buf_i]

    # Create output
dst, dst_lyr, dst_lyr_defn, field_names = Mk_outfile(out_shp.format(buffer_distance))
Mk_proj(utm33,out_shp.format(buffer_distance)[:-4])

for i in range(pts_lyr.GetFeatureCount()):
    pt_feature = pts_lyr.GetFeature(i)
    dst_lyr = Stamp_one_cookie(pt_feature, tile_idx_lyr, dst_lyr, dst_lyr_defn, field_names, buffer_distance)

dst = dst_lyr = None



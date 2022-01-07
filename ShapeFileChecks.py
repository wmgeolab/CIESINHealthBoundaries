import fiona
import shapely
from shapely.geometry import shape
from shapely.geometry import mapping, Polygon
import itertools
import geopandas as gpd
import csv
import os
import pandas as pd
import difflib
import unidecode
import topojson as tp
import matplotlib.pyplot as plt
import descartes
import pygeos

folder = "Intern_Completed/"
file = "Kasai_HealthZones"

reader = fiona.open(folder+file, 'r')

#reader.schema
#reader.driver
#reader.bounds

gdf = gpd.read_file(folder+file)
gdf.plot(figsize=(15,15), edgecolor="purple", facecolor="None")

#gdf.to_file("Converted_Files/Kasai_Oriental/Kasai_Oriental_Edited.shp")


error_log_name = "Error_Logs/" + file + "_Errors.txt"

with open(error_log_name, "w") as log:
  log.write("Errors:" + "\n")  

# Projection check
projection_error = None
projection = reader.crs
if projection['init'] == "epsg:4326":
  print("Projection Valid")
else:  
  if projection == {}:
    projection_error = "Projection Empty"
  else:
    projection_error = "Projection Invalid"
  with open(error_log_name, "a") as log:
    log.write(projection_error + "\n")



# copy set of shapes into new variable and store left and right property columns
new_shps = []
all_left = []
all_right = []
i = 0
for shp in reader:
  new_shps.append(shp)
  left = []
  right = []
  for l, r in shp['properties'].items(): 
    left.append(l)
    right.append(r)
  all_left.append(left)
  all_right.append(right) 
  i+=1


# Autofill zone & area ids

# define function to fill in shape values as blank when warranted
def fill_values(shp):
  try:
    shp['properties']['as_id']
  except KeyError:
    shp['properties']['as_id'] = ''
  try:
    shp['properties']['as_uid']
  except KeyError:
    shp['properties']['as_uid'] = ''
  try:
    shp['properties']['zs_id']
  except KeyError:
    shp['properties']['zs_id'] = ''
  try:
    shp['properties']['zs_uid']
  except KeyError:
    shp['properties']['zs_uid'] = ''



pyramid_fieldnames = ['dps_id','dps_uid','dps_province','zs_id','zs_uid',
                      'zone_de_sante','as_id','as_uid','aire_de_sante','fosa_id','fosa_uid']
pyramid = pd.read_csv(r"Copy of DSNIS_dhis2_pyramide_20200908 - Pyramide V2.csv", usecols=pyramid_fieldnames)


# for each shape, examine Aire de Sante and Zone de Sante for closest matches in pyramid
p=0
as_names = []
zs_names = []
for shp in new_shps:
  c=0
  # define AS and ZS indeces and names and values
  try:  
    AS_index = all_left[p].index('AS_')
    AS = all_right[p][AS_index]
    if not AS:
      raise ValueError
    as_names.append('AS_')
  except ValueError:
    try:
      AS_index = all_left[p].index('AS1_1')
      AS = all_right[p][AS_index]
      if not AS:
        raise ValueError
      as_names.append('AS1_1')
    except ValueError: 
      try:
        AS_index = all_left[p].index('aire_de_sa')
        AS = all_right[p][AS_index]
        if not AS:
          raise ValueError
      except ValueError:
        with open(error_log_name, "a") as log:
          log.write("Autofill: Aire de Sante not listed in source file. id: " + shp['id'] + "\n")
          fill_values(shp)
          as_names.append('N/A')  
          p+=1
          continue 
      as_names.append('aire_de_sa')
      
  try:
    ZS_index = all_left[p].index('ZS')
    ZS = all_right[p][ZS_index]
    if not ZS:
      raise ValueError
    zs_names.append('ZS')
  except ValueError:
    try:
      ZS_index = all_left[p].index('zone_de_sa')
      ZS = all_right[p][ZS_index]
      if not ZS:
        raise ValueError
    except ValueError:
      with open(error_log_name, "a") as log:
        log.write("Autofill: Zone de Sante not listed in source file. id: " + shp['id'] + "\n")
        fill_values(shp)
        zs_names.append('N/A')  
        p+=1
        continue 
    zs_names.append('zone_de_sa')
  
  AS = AS.lower()
  AS = unidecode.unidecode(AS)
  AS = AS.replace('-', ' ')
  AS_len = len(AS.split())
    
  ZS = ZS.lower()
  ZS = unidecode.unidecode(ZS)
  ZS = ZS.replace('-', ' ')
  ZS_len = len(ZS.split())

  all_row_match = []
  for index,row in pyramid.iterrows():
    ads_phrase = pyramid['aire_de_sante'][index]
    ads_phrase = ads_phrase.lower()
    ads_words = ads_phrase.split()
    if AS_len == 1:
      row_ratio = []
      for w in ads_words:
        ratio = difflib.SequenceMatcher(None,AS,w).ratio()
        row_ratio.append(ratio)
      ads_max_row_ratio = max(row_ratio)
    if AS_len == 2:
      segments = [' '.join([i,j]) for i,j in zip(ads_words, ads_words[1:])]
      row_ratio = []
      for w in segments:
        ratio = difflib.SequenceMatcher(None,AS,w).ratio()
        row_ratio.append(ratio)
      ads_max_row_ratio = max(row_ratio)
    if AS_len == 3:
      segments = [' '.join([i,j,k]) for i,j,k in zip(ads_words, ads_words[1:], ads_words[2:])]
      row_ratio = []
      for w in segments:
        ratio = difflib.SequenceMatcher(None,AS,w).ratio()
        row_ratio.append(ratio)
      ads_max_row_ratio = max(row_ratio)

    zds_phrase = pyramid['zone_de_sante'][index]
    zds_phrase = zds_phrase.lower()
    zds_words = zds_phrase.split()
    if ZS_len == 1:
      row_ratio = []
      for w in zds_words:
        ratio = difflib.SequenceMatcher(None,ZS,w).ratio()
        row_ratio.append(ratio)
      zds_max_row_ratio = max(row_ratio)
    if ZS_len == 2:
      segments = [' '.join([i,j]) for i,j in zip(zds_words, zds_words[1:])]
      row_ratio = []
      for w in segments:
        ratio = difflib.SequenceMatcher(None,ZS,w).ratio()
        row_ratio.append(ratio)
      zds_max_row_ratio = max(row_ratio)
    if ZS_len == 3:
      segments = [' '.join([i,j,k]) for i,j,k in zip(zds_words, zds_words[1:], zds_words[2:])]
      row_ratio = []
      for w in segments:
        ratio = difflib.SequenceMatcher(None,ZS,w).ratio()
        row_ratio.append(ratio)
      zds_max_row_ratio = max(row_ratio)
    
    row_match = (4/3)*ads_max_row_ratio+(2/3)*zds_max_row_ratio
    all_row_match.append(row_match)  
  
  row_max = max(all_row_match)
  
  if row_max<1.65:
    with open(error_log_name, "a") as log:
      log.write("Autofill: No close match in pyramid for Aire de Sante " + AS + " and Zone de Sante " + ZS + ". id: " + shp['id'] + "\n")
    fill_values(shp)
    p+=1
    continue
  possible_indeces = [x for x, y in enumerate(all_row_match) if y == row_max]
  if len(possible_indeces) == 1:
    pass
  else:
    u = 0
    mult = []
    for indx in possible_indeces:
      mult.append(pyramid['aire_de_sante'][indx])
      if u>0:
        if mult[u] == mult[u-1]:
          pass
        else:
          with open(error_log_name, "a") as log:
            log.write("Autofill: No close match in pyramid for Aire de Sante " + AS + " and Zone de Sante " + ZS + ". id: " + shp['id'] + "\n")
          c=1
          fill_values(shp)
          break
      u+=1
    if c == 1:
      p+=1
      continue
   
  py_index = possible_indeces[0]  
   
  AS_id = pyramid['as_id'][py_index]
  AS_uid = pyramid['as_uid'][py_index]
  ZS_id = pyramid['zs_id'][py_index]
  ZS_uid = pyramid['zs_uid'][py_index]
  
  try:
    shp['properties']['as_id']
  except KeyError:
    shp['properties']['as_id'] = AS_id
  else:
    if not shp['properties']['as_id']:
      shp['properties']['as_id'] = AS_id
  
  try:
    shp['properties']['as_uid']
  except KeyError:
    shp['properties']['as_uid'] = AS_uid
  else:
    if not shp['properties']['as_uid']:
      shp['properties']['as_uid'] = AS_uid
  
  try:
    shp['properties']['zs_id']
  except KeyError:
    shp['properties']['zs_id'] = str(ZS_id)
  else:
    if not shp['properties']['zs_id']:
      shp['properties']['zs_id'] = str(ZS_id)
  
  try:
    shp['properties']['zs_uid']
  except KeyError:
    shp['properties']['zs_uid'] = ZS_uid
  else:
    if not shp['properties']['zs_uid']:
      shp['properties']['zs_uid'] = ZS_uid

  p+=1




# Find invalid individual shapes

badgeoms = []
badgeomids = []
a=b=c=d=0
for shp in new_shps:
  geom = shape(shp['geometry'])
  if(not geom.is_valid):
    a += 1
    badgeoms.append(shp)
    if as_names[d] == 'N/A':
      a_s = 'N/A'
    else:
      a_s = gdf[as_names[d]][d]
    with open(error_log_name, "a") as log:
      log.write("Invalid individual polygon geometry for " + a_s + "\n")
    if(geom.buffer(0).is_valid):
      b += 1
    else:
      c += 1
  d+=1     

      
#    
# Mesh Checks
#


# overlap check
overlaps = []    
for index,row in gdf.loc[0:len(gdf)-2].iterrows():
  geom = shape(gdf['geometry'][index])
  for index2,row2 in gdf.loc[index+1:len(gdf)-1].iterrows():
    diff = geom.difference(shape(gdf['geometry'][index2]))
    if diff != geom:
      if as_names[index] == 'N/A':
        as1 = 'N/A'
      else:
        as1 = gdf[as_names[index]][index]
      if as_names[index2] == 'N/A':
        as2 = 'N/A'
      else:
        as2 = gdf[as_names[index2]][index2]
      with open(error_log_name, "a") as log:
        log.write("Overlap between " + as1 + " and " + as2 + "\n")

        



# gap check

bounds = reader.bounds
test_point = [(bounds[0],bounds[1])] 
distance_lon = bounds[2]-bounds[0]
distance_lat = bounds[3]-bounds[1]
bkgd = [shapely.geometry.Polygon([(x, y), (x + distance_lon, y), (x + distance_lon, y + distance_lat),(x, y + distance_lat)]) for x,y in test_point]
bkgd = bkgd[0]

for index,row in gdf.iterrows():
  geom = shape(gdf['geometry'][index])
  bkgd -= geom


gap_shape_path = 'Gaps_Shapes/'+ file + '_Gaps.shp'
gap_schema = {'geometry': '3D Polygon', 'properties': {'id':'int:10'}}
with fiona.open(gap_shape_path, 'w', 'ESRI Shapefile', gap_schema) as c:
    # ff there are multiple geometries, put the "for" loop here
    c.write({'geometry': mapping(bkgd), 'properties': {'id': 0}})
  

gap_image_path = 'Gaps_Images/'+ file + '_Gaps.png'   
gaps = gpd.read_file(gap_shape_path)
fig = gaps.plot(figsize=(15,15), edgecolor="blue", facecolor="black")     
fig = fig.get_figure()
fig.savefig(gap_image_path)





# create new shapefile

schema = reader.schema.copy()
try:
  schema['properties']['as_id']
except KeyError as e:
  schema['properties']['as_id']='str:100'
try:
  schema['properties']['as_uid']
except KeyError as e:
  schema['properties']['as_uid']='str:100'
try:
  schema['properties']['zs_id']
except KeyError as e:
  schema['properties']['zs_id']='str:100'
try:
  schema['properties']['zs_uid']
except KeyError as e:
  schema['properties']['zs_uid']='str:100'
  
    
final_filename = "Final_Files/geoBoundaries_" + file.replace(".shp","").replace("Intern_Completed/","")

with fiona.open(final_filename, 'w', reader.driver, schema, reader.crs) as output:
  for elem in new_shps:
    output.write(elem)



    
    

    
    
"""
# snapping
gdfTopo = tp.Topology(gdf, topology=True, 
                           #prequantize=False, 
                           shared_coords=True,
                           prevent_oversimplify=True).toposimplify(0.01).to_gdf().set_crs("EPSG:4326")
#gdfTopo.plot(figsize=(15,15)).get_figure().savefig("fixedMesh.png")

contrast = gpd.overlay(bkgd, gdf, how="symmetric_difference")
contrast.plot(figsize=(15,15)).get_figure().savefig("meshDifference.png")

#contrast.to_file("difference.geojson", driver='GeoJSON')

"""    
    
    
    
"""

The rules that define valid individual Polygons are:

a) Polygons are topologically closed;
b) The boundary of a Polygon consists of a set of LinearRings that make up its exterior and interior boundaries;
c) No two Rings in the boundary cross and the Rings in the boundary of a Polygon may intersect at a Point but
only as a tangent;
d) A Polygon may not have cut lines, spikes, or punctures;
e) The interior of every Polygon is a connected point set;
f) The exterior of a Polygon with 1 or more holes is not connected. Each hole defines a connected component of
the exterior.


"""


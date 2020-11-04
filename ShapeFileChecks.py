import fiona
from shapely.geometry import shape
import itertools
import geopandas as gpd
import csv
import os
import pandas as pd
import difflib
import unidecode
import itertools


file = "Intern_Completed/Lomami_Merged_Edits"

reader = fiona.open(file, 'r')

map1 = gpd.read_file(file)
map1.plot(figsize=(5,5), edgecolor="purple", facecolor="None")

#reader.schema
#reader.driver
#first = reader.next()
#first
#geom = shape(first['geometry'])
#geom
#geom.area
#list(geom.exterior.coords)
#list(geom.interiors)
#geom.geom_type
#len(geom.geoms)
#reader.bounds



error_log_name = "Error_Logs/" + file.replace(".shp","").replace("Intern_Completed/","") + " Errors.txt"

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
    
    
badgeoms = []
badgeomids = []
a=b=c=d=0
for feat in reader:
  geom = shape(feat['geometry'])
  a += 1
  if(not geom.is_valid):
    b += 1
    badgeoms.append(feat)
    with open(error_log_name, "a") as log:
      log.write("Invalid individual polygon geometry at ID: " + feat['id'] + "\n")
    if(geom.buffer(0).is_valid):
      c += 1
      #geom.difference(geom.buffer(0))
    else:
      d += 1




# Autofill zone & area ids

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
  
pyramid_fieldnames = ['dps_id','dps_uid','dps_province','zs_id','zs_uid',
                      'zone_de_sante','as_id','as_uid','aire_de_sante','fosa_id','fosa_uid']
pyramid = pd.read_csv(r"Copy of DSNIS_dhis2_pyramide_20200908 - Pyramide V2.csv", usecols=pyramid_fieldnames)




# for each shape, examine Aire de Sante and Zone de Sante for closest matches in pyramid
p=0
for shp in new_shps:
  c=0
  AS_index = all_left[p].index('AS_')
  AS = all_right[p][AS_index]
  if not AS:
    AS_index = all_left[p].index('AS1_1')
    AS = all_right[p][AS_index]
  AS = AS.lower()
  AS = unidecode.unidecode(AS)
  AS = AS.replace('-', ' ')
  AS_len = len(AS.split())
  ZS_index = all_left[p].index('ZS')
  ZS = all_right[p][ZS_index]
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
      log.write("Autofill: No close match in pyramid for Aire de Sante " + AS + " and Zone de Sante " + ZS + ". id: " + str(p) + "\n")
    try:
      shp['properties']['as_id']
    except KeyError as e:
      shp['properties']['as_id'] = ''
    try:
      shp['properties']['as_uid']
    except KeyError as e:
      shp['properties']['as_uid'] = ''
    try:
      shp['properties']['zs_id']
    except KeyError as e:
      shp['properties']['zs_id'] = ''
    try:
      shp['properties']['zs_uid']
    except KeyError as e:
      shp['properties']['zs_uid'] = ''
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
            log.write("Autofill: No close match in pyramid for Aire de Sante " + AS + " and Zone de Sante " + ZS + ". id: " + str(p) + "\n")
          c=1
          try:
            shp['properties']['as_id']
          except KeyError as e:
            shp['properties']['as_id'] = ''
          try:
            shp['properties']['as_uid']
          except KeyError as e:
            shp['properties']['as_uid'] = ''
          try:
            shp['properties']['zs_id']
          except KeyError as e:
            shp['properties']['zs_id'] = ''
          try:
            shp['properties']['zs_uid']
          except KeyError as e:
            shp['properties']['zs_uid'] = ''
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
  except KeyError as e:
    shp['properties']['as_id'] = AS_id
  else:
    if not shp['properties']['as_id']:
      shp['properties']['as_id'] = AS_id
  
  try:
    shp['properties']['as_uid']
  except KeyError as e:
    shp['properties']['as_uid'] = AS_uid
  else:
    if not shp['properties']['as_uid']:
      shp['properties']['as_uid'] = AS_uid
  
  try:
    shp['properties']['zs_id']
  except KeyError as e:
    shp['properties']['zs_id'] = str(ZS_id)
  else:
    if not shp['properties']['zs_id']:
      shp['properties']['zs_id'] = str(ZS_id)
  
  try:
    shp['properties']['zs_uid']
  except KeyError as e:
    shp['properties']['zs_uid'] = ZS_uid
  else:
    if not shp['properties']['zs_uid']:
      shp['properties']['zs_uid'] = ZS_uid

  p+=1


  

  
# create new shape file
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
  

  

final_filename = "Final_Shape_Files/geoBoundaries_" + file.replace(".shp","").replace("Intern_Completed/","")

with fiona.open(final_filename, 'w', 'ESRI Shapefile', schema, reader.crs) as output:
  for elem in new_shps:
    output.write(elem)



    
    
    
    
    
    
    



"""

The rules that define valid Polygons are:

a) Polygons are topologically closed;
b) The boundary of a Polygon consists of a set of LinearRings that make up its exterior and interior boundaries;
c) No two Rings in the boundary cross and the Rings in the boundary of a Polygon may intersect at a Point but
only as a tangent;
d) A Polygon may not have cut lines, spikes or punctures;
e) The interior of every Polygon is a connected point set;
f) The exterior of a Polygon with 1 or more holes is not connected. Each hole defines a connected component of
the exterior.


"""


#!pip3 install fiona
#!pip3 install geopandas
#!pip3 install geopandas descartes

import fiona
from shapely.geometry import shape
import itertools
import geopandas as gpd 

reader = fiona.open('geoBoundaries-3_0_0-AUS-ADM2-shp', 'r')

projection = reader.crs
projection


print(reader.schema)
first = reader.next()
print(first)
shp_geom = shape(first['geometry'])
print(shp_geom)
print(type(shp_geom))


bbox = reader.bounds
bbox
 
  
  
a=b=c=0
for feat in reader:
    geom = shape(feat['geometry'])
    a += 1
    if(not geom.is_valid):
      b += 1
      print(b)
      if(geom.buffer(0).is_valid):
        c + 1
        #Shows the errors:
        geom.difference(geom.buffer(0))
a  
b
c
geom.geom_type





mapshaper_reader = fiona.open('mapshaperclean', 'r')

a=0 
for (feat1,feat_2) in zip(reader,mapshaper_reader):
  a += 1
  geom1 = shape(feat1['geometry'])
  geom2 = shape(feat2['geometry'])
  if(geom1.difference(geom2).is_empty):
    print(a)

    
    
    
    
    
    
 


"""
green = reader
blue = mapshaper_reader 
# test the function difference between green and blue shapefiles
[not shape(i['geometry']).difference(shape(j['geometry'])).is_empty for i,j in zip(list(green),list(blue))]
[False, False, False, True]
# control
for geom in [shape(i['geometry']).difference(shape(j['geometry'])) for i,j in zip(list(green),list(blue))]:
  print(geom)
# test the function difference between blue and green shapefiles
[not shape(i['geometry']).difference(shape(j['geometry'])).is_empty for i,j in zip(list(blue),list(green))]
[True, False, False, False]
# control
for geom in [shape(i['geometry']).difference(shape(j['geometry'])) for i,j in zip(list(blue),list(green))]:
  print(geom)
# thus you can write a resulting shapefile withe the differences
from shapely.geometry import mapping
schema = {'geometry': 'Polygon','properties': {'test': 'int'}}
with fiona.open('diff.shp','w','ESRI Shapefile', schema) as e:
  for geom  in [shape(i['geometry']).difference(shape(j['geometry'])) for i,j in zip(list(green),list(blue))]:
    if not geom.is_empty:
      e.write({'geometry':mapping(geom), 'properties':{'test':1}})
  for geom  in [shape(i['geometry']).difference(shape(j['geometry'])) for i,j in zip(list(blue),list(green))]:
    if not geom.is_empty:
      e.write({'geometry':mapping(geom), 'properties':{'test':2}})
      
"""                                                                                      
                                                                                        
                                                                                        
                                                                                        
                                                                                        
                                                             
                                                                                       
map1 = gpd.read_file("geoBoundaries-3_0_0-AUS-ADM2-shp")
map1.plot(figsize=(5,5), edgecolor="purple", facecolor="None")

map2 = gpd.read_file("mapshaperclean")
map2.plot(figsize=(5,5), edgecolor="green", facecolor="None")

                                                                                                                 
#map1.difference(map2)  
  
 

  
###
##### CIESIN rar
###
  
  
CSNreader = fiona.open('DRC_healthboundaries_20200720.gdb', 'r')

projection = CSNreader.crs
projection


print(CSNreader.schema)
first = CSNreader.next()
print(first)
shp_geom = shape(first['geometry'])
print(shp_geom)
print(type(shp_geom))


bbox = CSNreader.bounds
bbox
 
  
  
a=b=c=0
for feat in CSNreader:
    geom = shape(feat['geometry'])
    a += 1
    if(not geom.is_valid):
      b += 1
      print(b)
      if(geom.buffer(0).is_valid):
        c + 1
        #Shows the errors:
        geom.difference(geom.buffer(0))
a  
b
c
geom.geom_type


CSNmap = gpd.read_file("DRC_healthboundaries_20200720.gdb")
CSNmap.plot(figsize=(5,5), edgecolor="purple", facecolor="None")






from shapely.geometry import Polygon
coords = [(0, 0), (0, 2), (1, 1), (2, 2), (2, 0), (1, 1), (0, 0)]
bowtie = Polygon(coords)
bowtie
bowtie.is_valid

clean = bowtie.buffer(0)
clean.is_valid
clean
list(clean[0].exterior.coords)
list(clean[1].exterior.coords)

diff = bowtie.difference(clean)
diff

from shapely.geometry import Point
a = Point(1, 1).buffer(1).buffer(-1)
a
b = Point(2, 1).buffer(1).buffer(-1)
b
a.difference(b)
b.difference(a)
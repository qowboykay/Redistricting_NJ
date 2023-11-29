#You must install these libraries onto your system before running the code
#!pip install geojson
#!pip install shapely
#!pip install PyShp

import geojson
import pandas as pd
import numpy as np
import time
import csv
import pickle

#Setting the state (NJ)
state = 'nj'

#assignment will contain the precinct-district assignments from my map
assignment = pd.read_csv('precinct-assignments-congress-'+state+'.csv')
precinct_list=[]
for p in range(0,len(assignment)):
    precinct_list.append(assignment.iloc[p].GEOID20)
len(precinct_list)
#print(precinct_list)

#Checking geography of maps
import shapefile as shp
from shapely.geometry import Polygon,shape,MultiPolygon


shpfile = 'Map_Data/'+state+'_vtd_2020_bound/'+state+'_vtd_2020_bound.shp'
dbffile = 'Map_Data/'+state+'_vtd_2020_bound/'+state+'_vtd_2020_bound.dbf'
shxfile = 'Map_Data/'+state+'_vtd_2020_bound/'+state+'_vtd_2020_bound.shx'
shpfile = shp.Reader(shp=shpfile, shx=shxfile, dbf=dbffile)
print(shpfile)

field_names = []
print(shpfile.fields[1:])
for f in shpfile.fields[1:]:
    field_names.append((f[0]))
print(field_names)

precinct_boundaries={}
count=0
for sr in shpfile.iterShapeRecords():
    geom = sr.shape # get geo bit
    rec = sr.record # get db fields
    precinct_boundaries[rec[3]]=geom
    count=count+1
    
##print the coordinates of a precinct polygon
if state=='nj':
    print(Polygon(shape(precinct_boundaries['34003060003'])))
if state=='nh':
    print(Polygon(shape(precinct_boundaries['33007SARG01'])))
###Examples
if state == 'nj':
    a = Polygon(shape(precinct_boundaries['34003060003']))
    b = Polygon(shape(precinct_boundaries['34007043046']))
    c = Polygon(shape(precinct_boundaries['34007043047']))

if state == 'nh':
    a = Polygon(shape(precinct_boundaries['33007RAND01']))
    b = Polygon(shape(precinct_boundaries['33007JEFF01']))
    c = Polygon(shape(precinct_boundaries['33007CARR01']))

#print(shape(precinct_boundaries['34029102002']))
                  
#b = Polygon(precint_boundaries['34007043046'])
print(a.touches(b))
print(c.touches(b))
print(c.touches(a))

#Function to check if two precinct overlap (needs precinct_boundaries to be instanciated)
def is_contiguous_precinct(p1,p2,precinct_boundaries):
    #print(shape(precinct_boundaries[p1]).type)
    try:
        if shape(precinct_boundaries[p1]).type == 'Polygon':
            a = Polygon(shape(precinct_boundaries[p1]))
        elif shape(precinct_boundaries[p1]).type == 'MultiPolygon':
            a = MultiPolygon(shape(precinct_boundaries[p1]))
        else:
            return False
        if shape(precinct_boundaries[p2]).type == 'Polygon':
            b = Polygon(shape(precinct_boundaries[p2]))
        elif shape(precinct_boundaries[p2]).type == 'MultiPolygon':
            b = MultiPolygon(shape(precinct_boundaries[p2]))
        else:
            return False
        return(a.touches(b))
    except KeyError: return False
    
if state == 'nj':
    print(is_contiguous_precinct('34003060003','34007043046',precinct_boundaries))
    print(is_contiguous_precinct('34007043047','34007043046',precinct_boundaries))
    print(is_contiguous_precinct('34003060003','34003060002',precinct_boundaries))
    print(is_contiguous_precinct('34007043047','34003060003',precinct_boundaries))
if state == 'nh':
    print(is_contiguous_precinct('33007CARR01','33007JEFF01',precinct_boundaries))
    print(is_contiguous_precinct('33007CARR01','33007RAND01',precinct_boundaries))

    #Function to find contiguous districts (needs precinct_boundaries to be instanciated)
    ## This code will tell you how many precinct are contiguous to a given precinct
def contiguous_precincts(p1,precinct_list,precinct_boundaries):
    count=0
    neighbors=[]
    for p in range(0,len(precinct_list)):
        #print(pre_data.iloc[p].GEOID20)
        if(is_contiguous_precinct(p1,precinct_list[p],precinct_boundaries)):
            count+=1
            neighbors.append(precinct_list[p])
    return([count,neighbors])

if state=='nh':
    print(contiguous_precincts('33007CARR01',precinct_list,precinct_boundaries))
if state=='nj':
    print(contiguous_precincts('34041080003',precinct_list,precinct_boundaries))

    ## NOTE: on my laptop, the code takes 33sec for NH, for NJ it takes 1h50m
tic = time.perf_counter()
    
contiguous={}
for p in range(0,len(precinct_list)):
    contiguous[precinct_list[p]]=[]

##Code is using pickling to be able to restart in case of crash. Needs to create a pickle directory first
## in case restart is needed, update the range of values for i
    
for i in range(0,len(precinct_list)):
    neighbors = contiguous_precincts(precinct_list[i],precinct_list[i+1:len(precinct_list)],precinct_boundaries)[1]
    #print(neighbors)
    contiguous[precinct_list[i]].extend(neighbors)
    for n in range(0,len(neighbors)):
        contiguous[neighbors[n]].append(precinct_list[i])
    if (i%100)==0:
        #save temp data in case it crashed
        picklename = "pickle/"+state+"_runAt"+str(i)+".p"
        with open(picklename, 'wb') as f:
            pickle.dump(contiguous, f) 
        toc = time.perf_counter()
        print(f"Pickle at " +str(i) + f",Took {toc - tic:0.4f} seconds") 
    
        
#save final data 
picklename = "pickle/"+state+"_FinalRun.p"
pickle.dump(contiguous, open(picklename, "wb")) 
toc = time.perf_counter()
print(f"Final Pickle at " +str(i) + f",Took {toc - tic:0.4f} seconds") 
        
with open('Contiguity_'+state+'.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    for key, value in contiguous.items():
       writer.writerow([key, value])
toc = time.perf_counter()

print(f"Took {toc - tic:0.4f} seconds")

#Make sample maps to test the contiguity mappings

picklename = "pickle/"+state+"_FinalRun.p"
contiguous_precincts= pickle.load(open(picklename,"rb"))

if state=='nj':
    print(contiguous_precincts['34041080003'])
if state=='nh':
    print(contiguous_precincts['33007JEFF01'])

clear_map=pd.DataFrame(precinct_list)
clear_map[1]=np.repeat(1,len(precinct_list))

#create a random map by selecting n precincts and coloring their neighbors
import random

n=5

new_map=pd.DataFrame(precinct_list)
new_map[1]=np.repeat(1,len(precinct_list))
new_map.columns=['A','B']
index=2
sample = random.sample(precinct_list,n)
for s in range(0,len(sample)):
    print(sample[s])
    neighbors = contiguous_precincts[sample[s]]
    #print(new_map.loc[new_map['A'] == sample[s],])
    new_map.loc[new_map['A'] == sample[s],'B']=index
    for n in range(0,len(neighbors)):
        new_map.loc[new_map['A'] == neighbors[n],'B']=index+1
    index=index+2
new_map.to_csv("pickle/testmap1.csv",index=False)
                                         
        


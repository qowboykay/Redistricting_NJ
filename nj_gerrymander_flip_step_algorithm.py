#You must install these libraries onto your system before running the code
#!pip install geojson
#!pip install shapely
#!pip install PyShp
#!pip install networkx

import geojson
import pandas as pd
import numpy as np
import networkx as nx
import time
import csv
import ast
import shapefile as shp
from shapely.geometry import Polygon,shape,MultiPolygon
import shapely.ops
import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

#Helper functions

def isDistrictContiguous(district_num, assignment, contiguity_list, print_isolates=False, ignore_list=[]):
    ## input:
    ## district_num: the district number
    ## assignment: the assignment from precinct to district
    ## contiguity_list: the list of neighbors for each precinct, from the csv file
    contiguity_list.columns = ['Precinct','Neighbors']
    district_graph = nx.Graph() #creates an empty undirected graph
    district_nodes = assignment[assignment['District']==district_num]['GEOID20'].tolist()
    for i in ignore_list:
        try:
            district_nodes.remove(i)
        except ValueError:
            pass
    district_graph.add_nodes_from(district_nodes)
    for id in district_nodes:
        neighbors = ast.literal_eval(contiguity_list[contiguity_list['Precinct']==id]['Neighbors'].values.tolist()[0])
        # needed to convert string to list because the csv encodes the list as a string
        for neighbor in neighbors:
            district_graph.add_edge(id,neighbor)
    if(print_isolates):
        print(list(nx.isolates(district_graph)))
    return nx.is_connected(district_graph)

def getDistrictPopulations(assignment,data_file, num_district):
    population = {}
    for i in range (1,num_district+1):
        population[i] = data_file[data_file['GEOID20'].isin(assignment[assignment['District']==i]['GEOID20'])]['Total_2020_Total'].sum()
    return population

def getDistrictShape(district_id, assignment, boundaries):
    list_precincts = assignment[assignment['District']==district_id]['GEOID20']
    precinct_shapes = []
    for i in list_precincts:
        if shape(boundaries[i]).type == 'Polygon':
            precinct_shapes.append(Polygon(shape(boundaries[i])))
        elif shape(boundaries[i]).type == 'MultiPolygon':
            precinct_shapes.append(MultiPolygon(shape(boundaries[i])))
    district_shape = shapely.ops.unary_union(precinct_shapes)
    #print(district_shape)
    return district_shape

def pp_compactness(geom): # Polsby-Popper
    p = geom.length
    a = geom.area
    return (4*np.pi*a)/(p*p)

def box_reock_compactness(geom): # Reock on a rectangle bounding box
    a = geom.area
    bb = geom.bounds # bounds gives you the minimum bounding box (rectangle)
    bba = abs(bb[0]-bb[2])*abs(bb[1]-bb[3])
    return a/bba

# This is the current assignment of precinct to congressional districts (12 of them for NJ)"""

nj_current_assignment = pd.read_csv('/content/sample_data/precinct-assignments-congress-nj.csv')
nj_current_assignment

# This is the current demographic and voter data
#The data has a lot of attributes that lists voters of different demographics and parties in different elections.
#For this code I will only keep votes from the 2020  presidential election and the total 2020 population counts.
#You can use additional columns (e.g., Governor's elections results, voting age (VAP) population counts, or the composite Dem/Rep score)

nj_precinct_data = pd.read_csv('/content/sample_data/precinct-data-congress-nj.csv')
keepcolumns = ['GEOID20','District','Total_2020_Pres','Dem_2020_Pres','Rep_2020_Pres','Total_2020_Total','White_2020_Total','Hispanic_2020_Total','Black_2020_Total','Asian_2020_Total','Native_2020_Total','Pacific_2020_Total']
nj_precinct_data = nj_precinct_data[keepcolumns]
nj_precinct_data

# This is the precinct boundary data (uses shapely)
#This is data that represents the geography of the districts. It is needed to test for contiguity, or for any districting partitioning method based on geography.
#The data is in Shapely format. Each district is represented as a set of points that are connected to create the district shape (in the long/lat coordinates).
#Shapely geometric functions can be used to compare the shapes.
#These can be quite inefficient to run, so I pre-computed an index that, for each district, lists the districts that are contiguous to it. The code to generate the index is in Contiguity.py
#To manipulate the shapes, cast them into Shapely Polygons (see example below) and you can use the Polygon properties and functions:
#https://shapely.readthedocs.io/en/stable/reference/shapely.Polygon.html#shapely.Polygonshpfile = '/content/sample_data/nj_vtd_2020_bound.shp'

dbffile = '/content/sample_data/nj_vtd_2020_bound.dbf'
shxfile = '/content/sample_data/nj_vtd_2020_bound.shx'


shpfile = shp.Reader(shp=shpfile, shx=shxfile, dbf=dbffile)
nj_precinct_boundaries={}
for sr in shpfile.iterShapeRecords():
    geom = sr.shape # get geo bit
    rec = sr.record # get db fields
    nj_precinct_boundaries[rec[3]]=geom

# This is the precinct boundary data

#This uses the contiguity index I have pre-computed using Contiguity.py, that is stored in Contiguity_nj.csv.


nj_contiguity = pd.read_csv('/content/sample_data/Contiguity_nj.csv', header=None)

for i in range(1,13):
    print("District "+str(i)+" "+str(isDistrictContiguous(12, nj_current_assignment, nj_contiguity)))

#Compactness of the current assignment
for district in range(1,13):
    print("D"+str(district)+" PP : "+str(pp_compactness(getDistrictShape(district,nj_current_assignment,nj_precinct_boundaries))))
    print("D"+str(district)+" BR : "+str(box_reock_compactness(getDistrictShape(district,nj_current_assignment,nj_precinct_boundaries))))

# District Population of the current assignment
print(getDistrictPopulations(nj_current_assignment,nj_precinct_data, 12))


nj_flipstep_assignment = nj_current_assignment.copy()

# Function to calculate the political leaning of a district
def calculate_political_leaning(district_num, assignment, data_file):
    district_voters = data_file[data_file['GEOID20'].isin(assignment[assignment['District'] == district_num]['GEOID20'])]
    dem_voters = district_voters['Dem_2020_Pres'].sum()
    rep_voters = district_voters['Rep_2020_Pres'].sum()
    return dem_voters, rep_voters

# Function to flip border precincts with political and demographic considerations
def flip_border_precincts(district_from, district_to, max_white_percentage):
    border_precincts = []
    district_list = nj_flipstep_assignment[nj_flipstep_assignment['District'] == district_from]['GEOID20']

    for precinct in district_list:
        neighbors = ast.literal_eval(nj_contiguity[nj_contiguity['Precinct'] == precinct]['Neighbors'].values.tolist()[0])
        for neighbor in neighbors:
            if nj_flipstep_assignment.loc[nj_flipstep_assignment['GEOID20'] == neighbor, 'District'].values[0] == district_to:
                border_precincts.append(precinct)
                break

    # Randomly sample precincts and flip them unless it breaks contiguity or white population constraint
    flipped_precincts = np.random.choice(border_precincts, min(10, len(border_precincts)), replace=False)
    for flip in flipped_precincts:
        # Temporarily flip the precinct
        nj_flipstep_assignment.loc[nj_flipstep_assignment['GEOID20'] == flip, 'District'] = district_to

        # Check contiguity
        if not isDistrictContiguous(district_from, nj_flipstep_assignment, nj_contiguity) or \
           not isDistrictContiguous(district_to, nj_flipstep_assignment, nj_contiguity):
            nj_flipstep_assignment.loc[nj_flipstep_assignment['GEOID20'] == flip, 'District'] = district_from
            print("Contiguity broken by flipping precinct " + flip)
            continue

        # Check white population percentage
        district_population = getDistrictPopulations(nj_flipstep_assignment, nj_precinct_data, 12)[district_to]
        district_white_population = nj_precinct_data[nj_precinct_data['GEOID20'].isin(nj_flipstep_assignment[nj_flipstep_assignment['District'] == district_to]['GEOID20'])]['White_2020_Total'].sum()
        if (district_white_population / district_population) > max_white_percentage:
            nj_flipstep_assignment.loc[nj_flipstep_assignment['GEOID20'] == flip, 'District'] = district_from
            print("White population percentage exceeded by flipping precinct " + flip)

        # Additional check to favor Democrats
        dem_voters, rep_voters = calculate_political_leaning(district_to, nj_flipstep_assignment, nj_precinct_data)
        if dem_voters < rep_voters:
            nj_flipstep_assignment.loc[nj_flipstep_assignment['GEOID20'] == flip, 'District'] = district_from
            print("Flipping precinct " + flip + " does not favor Democrats. Reverting.")

# Maximum white population percentage, inversely affects minority representation
max_white_percentage = 0.72  

# Apply the flip-step algorithm iteratively over district pairs
for district_from in range(1, 13):
    for district_to in range(1, 13):
        if district_from != district_to:
            flip_border_precincts(district_from, district_to, max_white_percentage)

# Check contiguity for all districts after flipping
for district in range(1, 13):
    contiguity_status = isDistrictContiguous(district, nj_flipstep_assignment, nj_contiguity)
    print(f"District {district} contiguity: {contiguity_status}")


nj_flipstep_assignment.to_csv('nj_flipstep_map_v2.csv', index=False)

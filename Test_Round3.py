#Test_Round3.py



class coordinate(object):
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z


#Building has eight coordinates and a center coordinate and Width, Length and Height
class building(object):

    def __init__(self, eight_Coordinates,bottom_plane=[],top_plane=[]):
        self.eight_Coordinates = eight_Coordinates
        self.bottom_plane = []
        self.top_plane = []
        #print eight_Coordinates[0].x

        for coord in self.eight_Coordinates:
            if coord.z == 0:
                #print "bottom plane"
                self.bottom_plane.append(coord)
                #print self.bottom_plane
            else:
                self.top_plane.append(coord)

    #Task 1    
    def get_Coordinates(self):
        return self.eight_Coordinates	
    
    #Task 2
    def get_Center_of_Building(self):
        '''Break this problem down into two smaller problems
           1-first find the center on the bottom plane
           2-then find half the height of the building  
           then combine both those point to make the center'''

        ###The center of the bottom plane should be the mid-point of the diagnol of the bottom plane
        ###To find diagnol of the bottom plane just pick any of the four corners of bottom plane and calculate which of the other three points has the longest distance
        diagnol_length=0
        index_of_second_diagnol_point = 0 #the first index I chose as zero
        for i in range(1,4):
            
            tmplength = ((self.bottom_plane[0].x - self.bottom_plane[i].x)**2 + (self.bottom_plane[0].y - self.bottom_plane[i].y)**2)**0.5
            
            if tmplength > diagnol_length:
                diagnol_length = tmplength
                index_of_second_diagnol_point = i
        '''
        summary until now-->
        The list bottom_plane holds 4 coordinates of the bottom plane 
        I chose the 1st coordinate in that list to be the first point of the diagnol
        And the the second coordinate of that diagnol is the (index_of_second_diagnol_point)th coordinate in that list
        '''
        
        x_c = (abs(self.bottom_plane[0].x - self.bottom_plane[index_of_second_diagnol_point].x)/2) + min(self.bottom_plane[0].x,self.bottom_plane[index_of_second_diagnol_point].x)
        y_c = (abs(self.bottom_plane[0].y - self.bottom_plane[index_of_second_diagnol_point].y)/2) + min(self.bottom_plane[0].y,self.bottom_plane[index_of_second_diagnol_point].y)
        #Find half the height
        z_c = self.top_plane[0].z / 2
        
        #create the coordinate 
        center_coordinate = coordinate(x_c,y_c,z_c)
        
        #Return the full coordinate
        return center_coordinate
    

    #Task 3

    def get_Width(self):
    	diagnol_length=0
        index_of_second_diagnol_point = 0 
        for i in range(1,4):
            tmplength = ((self.bottom_plane[0].x - self.bottom_plane[i].x)**2 + (self.bottom_plane[0].y - self.bottom_plane[i].y)**2)**0.5
            if tmplength > diagnol_length:
                diagnol_length = tmplength
                index_of_second_diagnol_point = i
        indexes_of_adjacent_coords = [1,2,3]
        indexes_of_adjacent_coords.remove(index_of_second_diagnol_point) 

        # Now indexes_of_adjacent_coords holds the two coordinates of the two points on the bottom plane adjacent to the coordinate in index zero

        width = ((self.bottom_plane[0].x - self.bottom_plane[indexes_of_adjacent_coords[0]].x)**2 + (self.bottom_plane[0].y - self.bottom_plane[indexes_of_adjacent_coords[0]].y)**2)**0.5

        return width


    def get_Length(self):

    	diagnol_length=0
        index_of_second_diagnol_point = 0 
        for i in range(1,4):
            tmplength = ((self.bottom_plane[0].x - self.bottom_plane[i].x)**2 + (self.bottom_plane[0].y - self.bottom_plane[i].y)**2)**0.5
            if tmplength > diagnol_length:
                diagnol_length = tmplength
                index_of_second_diagnol_point = i
        indexes_of_adjacent_coords = [1,2,3]
        indexes_of_adjacent_coords.remove(index_of_second_diagnol_point) 

        # Now indexes_of_adjacent_coords holds the two coordinates of the two points on the bottom plane adjacent to the coordinate in index zero

        length = ((self.bottom_plane[0].x - self.bottom_plane[indexes_of_adjacent_coords[1]].x)**2 + (self.bottom_plane[0].y - self.bottom_plane[indexes_of_adjacent_coords[1]].y)**2)**0.5

        return length

    def get_Height(self):

        return self.top_plane[0].z - self.bottom_plane[0].z



#import library for parsing xml files and call it ET
import xml.etree.ElementTree as ET

#All of the input file is extracted and parsed and stored as xml objects in "tree"
tree = ET.parse("Test_Round3.dae")

#Input/output file
fileIn = open("Test_Round3.dae")
fileOut = open("Test_Round3_Kozuch.txt",'w')

#Task 1------------------------------------------------------------------------------
wholefile = fileIn.read()
numBuildings = wholefile.count("<geometry")
fileOut.write("Task 1:\n" + "There are "+str(numBuildings)+" Buildings in this .dae file\n\n")



#Task 2------------------------------------------------------------------------------
#The root in this case in the highest level TAG called Collada
root = tree.getroot()

#List to store all buildings
buildings = []

for i in range(0,numBuildings):
    #Extract the coordinates
    list_of_72Numbers = root[3][i][0][0][0].text.split(" ")

    #Get Unique Coordinates for each set of 72 numbers and Put them in vector format
    list_of_coordinates = []
    for i in range(0,72,3):
        coord = coordinate(float(list_of_72Numbers[i])*0.0254,float(list_of_72Numbers[i+1])*0.0254,float(list_of_72Numbers[i+2])*0.0254)
    
        unique = True                   #Assume its unique
        for tmp in list_of_coordinates: #loop through the existing unique ones
            if (tmp.x == coord.x) and (tmp.y == coord.y) and (tmp.z == coord.z):
                unique = False          #this coord is not unique

        if unique:
            list_of_coordinates.append(coord)

    #print list_of_coordinates

    
    buildings.append(building(list_of_coordinates))

    

#Task 2
###############################################################################################
fileOut.write("All coordinates in form (x,y,z) where z is up and x,y are the horizontal plane\n\n")
fileOut.write("Task 2:\n")

for i in range(0,numBuildings):
    fileOut.write("\nBuilding number "+ str(i+1)+"'s coordinates\n")
    eight_Coordinates= buildings[i].get_Coordinates()
    for coord in eight_Coordinates:
        fileOut.write("("+str(coord.x)+", "+str(coord.y)+", "+str(coord.z)+")\n")



#Task 3
###############################################################################################
fileOut.write("\n\nTask 3:\n")

for i in range(0,numBuildings):
    #print i	
    fileOut.write("\nBuilding number "+ str(i+1)+"'s center coordinate\n")
    center_coord = buildings[i].get_Center_of_Building()
    #print center_coord
    fileOut.write("("+str(center_coord.x)+", "+str(center_coord.y)+", "+str(center_coord.z)+")\n")

    
#Task 4
###############################################################################################
fileOut.write("\n\nTask 4:\n")

for i in range(0,numBuildings):
    fileOut.write("\nBuilding number "+ str(i+1)+"'s Length, width and height\n")
    fileOut.write("width --> " + str(buildings[i].get_Width()))
    fileOut.write("\nlength--> " + str(buildings[i].get_Length()))
    fileOut.write("\nheight--> " + str(buildings[i].get_Height()))


#Task 5
###############################################################################################
fileOut.write("\n\nTask 5:\n")

for i in range(0,numBuildings):
    w = buildings[i].get_Width()
    l = buildings[i].get_Length() 
    h = buildings[i].get_Height()

    surfaceArea = 2*w*l 
    surfaceArea+= 2*h*w
    surfaceArea+= 2*h*l 
    fileOut.write("\nBuilding number "+ str(i+1)+"'s surface area\n")
    fileOut.write("Surface area--> " + str(surfaceArea))


#Task 6
###############################################################################################
fileOut.write("\n\nTask 6:\n")

for i in range(0,numBuildings):
    w = buildings[i].get_Width()
    l = buildings[i].get_Length() 
    h = buildings[i].get_Height()

    volume = l*w*h
    fileOut.write("\nBuilding number "+ str(i+1)+"'s volume\n")
    fileOut.write("volume--> " + str(volume))


#Task 7
###############################################################################################

import math
def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)


fileOut.write("\n\nTask 7:\n")
'''
for i in range(0,nCr(numBuildings,2)):
    for j in range(0,nCr(numBuildings,2)).remove(i):
        print j
'''
print range(0,5).type












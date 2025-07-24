import sys
import pandas
import shapely
import numpy
import random

d = pandas.read_csv('popComp.csv' ,  on_bad_lines='skip' , delimiter=';')
d = d.values.tolist()

popRef = []

for i in range(len(d)):
    popRef_i = {}
    popRef_i['ID'] = d[i][0]
    popRef_i['geometry'] = d[i][1]
    popRef.append(popRef_i)


def enleveVirguleZ(string):
    char = ''
    chaine = ''
    for i in range(len(string)):
        
        if string[i] == ' ':
            chaine += char 
            char =' '
        
        elif string[i] == ',':
            char =''
        
        else : char += string[i]
        
    return chaine


def readPolygon(polygon):
    
    listPoint = polygon[12:len(str(polygon))]
    
    for i in range(len(listPoint)):
        if listPoint[i] == "(":
            return []
    
    listPoint = enleveVirguleZ(listPoint)
        
    
    
    lat , lon = '' , ''
    coord = 'lat'
    Point = []
    
    for i in range(len(listPoint)):
        #lat 
        if coord == 'lat' : 
            lat += listPoint[i]
            if listPoint[i] == " ": 
                coord = 'lon'
        else : 
            lon += listPoint[i]
            if listPoint[i] == " ": 
                coord = 'lat'
                Point.append((float(lat), float(lon)))
                lat = ''
                lon = ''
    
    Point.append(Point[0])
    
    return Point

def centreDesMassesP(polygon):
    
    x = 0
    y = 0
    
    for i in range(len(polygon)):
        x += polygon[i][0]
        y += polygon[i][1]
        
    return ( x/len(polygon) , y/len(polygon) )
    
    

def first(pop , a):
    popComp = []
    for i in range(len(pop)):
        proba = random.random()
        decallage = random.random()
        popComp_i = {}
        popComp_i['ID'] = pop[i]['ID']
        Points = pop[i]['geometry']
        if proba < a  : 
            geometry = []
            for i in range(len(Points)):
                geometry.append((Points[i][0] + decallage , Points[i][1] + decallage))
            
            if len(geometry) == 0 : 
                popComp_i['geometry'] = pop[i]['geometry']
                
            else : popComp_i['geometry'] = geometry

        else : popComp_i['geometry'] = pop[i]['geometry']
        popComp.append(popComp_i)
        
        
    return popComp 
    

    
def second(listGeometry, pop , bruiteur , echelle , spliteur ):
    
    listGeometry2 = []
    popComp = []    
    
    
    
    for i in range(len(pop)):
    
        proba = random.random()
        
        popComp_i = {}
        Points = pop[i]['geometry']
        
        if proba < bruiteur :
            
            
            #bruiteur
            geometry = []
            for j in range(len(Points)):
                distance = random.random()
                proba_decal = random.random()
                if  proba_decal < 0.25 : 
                    geometry.append((Points[j][0] + distance , Points[j][1] + distance))
                elif proba_decal > 0.25 and proba_decal < 0.5 : 
                    geometry.append((Points[j][0] + distance , Points[j][1] - distance))
                elif proba_decal > 0.5 and proba_decal < 0.75 : 
                    geometry.append((Points[j][0] - distance , Points[j][1] + distance))
                else :
                    geometry.append((Points[j][0] - distance , Points[j][1] - distance))
                
            popComp_i['ID'] = pop[i]['ID']
            popComp_i['geometry'] = geometry
            popComp.append(popComp_i)
            
            listGeometry2.append((listGeometry[i] , geometry ))
            
            
            
        if proba > bruiteur and proba < echelle : 
            #echelle 
            geometry = []
            c = centreDesMassesP(Points)
            dist = 0.5
            for j in range(len(Points)):
                PointA = Points[j]
                PointB = c
                
                if PointA[0] - PointB[0] == 0 : a = 0
                else : a = (PointA[1] - PointB[1]) / (PointA[0] - PointB[0])
                b = PointA[1] - a*PointA[0]
                
                dAB = ((PointA[0] - PointB[0])**2  + (PointA[1] - PointB[1])**2 )**.5
                
                if dAB == 0 : 
                    popComp_i['ID'] = pop[i]['ID']
                    popComp_i['geometry'] = pop[i]['geometry']
                    popComp.append(popComp_i)
                    
                else :
                    
                    dBC = abs(PointA[0] - PointB[0])
                    theta = numpy.arccos(dBC / dAB)
                    
                    x = PointA[0]
                    y = PointA[1]
                    
                    if x < c[0] and y > c[1] : 
                        
                        X = x + numpy.cos(theta)*dist
                        Y = a * X + b
                        
                    elif x < c[0] and y < c[1] : 
                        
                        X = x + numpy.cos(theta)*dist
                        Y = a * X + b
                        
                    elif x > c[0] and y > c[1] : 
                        
                        X = x - numpy.cos(theta)*dist
                        Y = a * X + b
                        
                    else :
                        
                        X = x - numpy.cos(theta)*dist
                        Y = a * X + b
                    
                    geometry.append((X,Y))
            
            popComp_i['ID'] = pop[i]['ID']
            popComp_i['geometry'] = geometry
            popComp.append(popComp_i)
            
            listGeometry2.append((listGeometry[i] , geometry ))
        
        #split
        if proba > echelle and proba < spliteur : 
            
            centroid = centreDesMassesP(Points)
            
            a = random.randrange(0,2)
            b = random.random()
            if b < 0.5 : 
                X = centroid[0] + a
                Y = centroid[1] + a
            else : 
                X = centroid[0] - a
                Y = centroid[1] - a
            
            Poi = []
            for j in range(len(Points)):
                Poi.append(Points[j])
            Poi.append(Points[1])
            Poi.append(Points[2])
    
            geometry = []
    
            for j in range(len(Points)):
                
                PointA = Poi[j]
                PointB = Poi[j+1]
                PointC = Poi[j+2]
                
                if PointA[0] - PointB[0] == 0 : a = 0
                else : a = (PointA[1] - PointB[1]) / (PointA[0] - PointB[0])
                b = PointA[1] - a*PointA[0]
                
                x = (X + a*(Y - b))/(1 + a**2)
                y = a*x + b
                
                pas = 0
                
                if PointA[1] < PointB[1] : 
                    if PointA[0] < x and x < PointB[0] : 
                        pas = 1
                else : 
                    if PointB[0] < x and x < PointA[0] : 
                        pas = 1
                
                if PointA[1] < PointB[1] : 
                    if PointA[1] < y and y < PointB[1] : 
                        pas = 1
                else : 
                    if PointB[1] < y and y < PointA[1] : 
                        pas = 1
                        
                if pas == 0 : 
                    geometry.append([])
                
                if PointB[0] - PointC[0] == 0 : aC = 0
                else : aC = (PointB[1] - PointC[1])/ (PointB[0]- PointC[0])
                bC = PointB[1] - aC*PointB[0]
                
                xC = (X + aC*(Y - bC))/(1 + aC**2)
                yC = aC*xC + bC
                
                pas = 0
                
                if PointC[1] < PointB[1] : 
                    if PointC[0] < xC and xC < PointB[0] : 
                        pas = 1
                else : 
                    if PointB[0] < xC and xC < PointC[0] : 
                        pas = 1
                
                if PointC[1] < PointB[1] : 
                    if PointC[1] < yC and yC < PointB[1] : 
                        pas = 1
                else : 
                    if PointB[1] < yC and yC < PointC[1] : 
                        pas = 1
                        
                if pas == 0 : 
                    geometry.append([])
                
                geometry.append(( (x,y) , PointB , (xC , yC) , (X,Y) , (x,y)))
                
            for j in range(len(geometry)):
                bati = {}
                bati['ID'] = str(pop[i]['ID']) + str(j)
                bati['geometry'] = geometry[j]
                popComp.append(bati)
                
                listGeometry2.append((listGeometry[i] , geometry[j] ))
            
        if proba > spliteur : 
            
            popComp_i['ID'] = pop[i]['ID']
            popComp_i['geometry'] = pop[i]['geometry']
            popComp.append(popComp_i)
            
            listGeometry2.append((listGeometry[i] , listGeometry[i] ))
        
    #return popComp
    return listGeometry2
    
def ter(popRef, popComp , a ):
    
    popRefN = []
    popCompN = []
    
    for i in range(len(popRef)):
        if len(popRef[i]) == 0 : 
            continue 
        
         
        if random.random() < 1 - a :
            popRefN.append(popRef[i])
            
    for i in range(len(popComp)):
        
        if len(popComp[i]) == 0 : 
            continue 
        
        if random.random() < 1 - a :
            popCompN.append(popComp[i])
        
            
    return (popRefN, popCompN)
    
def ter4geometry(listGeometry2 , a ):

    listGeometry2Final = []
    
    for i in range(len(listGeometry2)):
        if len(listGeometry2[i]) == 0 : 
            continue 
        
        if random.random() < 1 - a :
            listGeometry2Final.append(listGeometry2[i])
            
        else : 
            if random.random() > 0.5 :
                listGeometry2Final.append((listGeometry2[i][0] , [] ))
            else : 
                listGeometry2Final.append(([] , listGeometry2[i][1] ))
        
    return listGeometry2Final

    
def main(pop):
    
    popComp = [] 
    for i in range(len(pop)):
        popComp_i = {}
        popComp_i['ID'] = pop[i]['ID']
        popComp_i['geometry'] = readPolygon(pop[i]['geometry'])
        if len(popComp_i['geometry']) != 0:
            popComp.append(popComp_i)
    
    a = decallage
    
    firstL = first(popComp , a)
    
    b = bruiteur
    c = echelle
    d = spliteur
    
    if b > c or c > d or b > d :
        density = []
        
    else : 
        listGeometry = [ readPolygon(pop[i]['geometry']) for i in range(len(pop)) ]
    
        #secondL = second(listGeometry, firstL , b , c , d )
        
        listeFinal = second(listGeometry, firstL , b , c , d )
        
        e = suppression
        
        listeFinal = ter4geometry(listeFinal, e )
        
        #calcul distance surfacique 
        dSurfacique = []
        for i in range(len(listeFinal)):
            
            if len(listeFinal[i][0]) == 0 or len(listeFinal[i][1]) == 0 :
                # a etudier 
                dSurfacique.append(0)
            
            if len(listeFinal[i][0]) < 3 or len(listeFinal[i][1]) <3 : 
                continue
                
            geomRef = shapely.Polygon(listeFinal[i][0]).buffer(0)
            geomComp = shapely.Polygon(listeFinal[i][1]).buffer(0)
            
            if geomComp.is_valid is False  :
                continue
            else : 
                inter = shapely.intersection( geomRef , geomComp)
                union = shapely.union(geomRef,geomComp)
                #if union.area == 0 : 
                #    dSurfacique.append(1) 
                #else : 
                dSurfacique.append( 1 - inter.area /union.area )
        
        
        density = []
        for j in range(20):
            s = 0
            for i in range(len(dSurfacique)):
                if dSurfacique[i] >= 0.05*j and dSurfacique[i] < 0.05*(j+1) : 
                    s+=1
            density.append(s)
            
        s = 0
        for i in range(len(dSurfacique)): 
            if dSurfacique[i] == 1 : 
                s += 1
                
        density[19] += s
        
        for i in range(len(density)) : 
            density[i] /= len(dSurfacique)

    return density

density = main(popRef)

def calcul(listeA , listeB) :

    if len(density) == 0 :
    
        return 100
        
    sum_ = 0
    
    for i in range(len(listeA)):
        
        sum_ += abs( listeA[i] - listeB[i]) 
        
    return sum_


listeB = [0.3594998263285863, 0.013546370267453978, 0.031087183049670026, 0.03143452587704064, 0.007467870788468218, 0.007120527961097603, 0.004862799583188607, 0.004168113928447377, 0.0019103855505383813, 0.0005210142410559222, 0.0001736714136853074, 0.0001736714136853074, 0.0015630427231677665, 0.0006946856547412296, 0.013546370267453978, 0.016151441472733587, 0.0006946856547412296, 0.0017367141368530739, 0.0013893713094824591, 0.5022577283779089]

calcul = calcul(density , listeB)

import math

if calcul != 100 : 
    diff = calcul
else : 
    diff = math.nan

print('Fin de python task')

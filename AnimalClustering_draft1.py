#!/usr/bin/python
#Christian Felt
#Data Mining, spring 2017
#Performs clustering in various ways
#on spatial data from the IUCN Red List of Endangered Species. 
#The data must be read into matrices where the first column
#is latitude and the second column is longitude coordinates
#representing the centroid of each polygon in the original
#spatial data. (I used ArcMap software as well as the code
#here to preprocess the data.) 

import numpy as np
import matplotlib.pyplot as plt
import random
import copy
import math
import csv
from math import radians, cos, sin, sqrt, atan2
#from matplotlib import mlab

#In my csv files, the lat and long are in columns 27 and 28
def readCSVFileIntoMatrix(filename, latitudeColumn, longitudeColumn):
	
	#First we count the number of rows in the csv file
	numRows = 0
	with open(filename, 'r') as csvFile:
		dataReader = csv.reader(csvFile, delimiter=',')
		for row in dataReader:
			numRows += 1
	
	dataMatrix = np.zeros([numRows, 2])

	with open(filename, 'r') as csvFile:
			dataReader = csv.reader(csvFile, delimiter=',')
			i = 0
			for row in dataReader:
				if (i==0):
					i += 1
				else:
					longitude = row[27]
					latitude = row[26]
					dataMatrix[i][0] = latitude
					dataMatrix[i][1] = longitude
					i += 1
				
	return dataMatrix	

#Computes the great-circle distance between two points
#on the Earth's surface using the Haversine formula as described here:
#http://www.movable-type.co.uk/scripts/latlong.html
#We model the Earth as a sphere, which can lead to errors up to 0.3% or 22km
def greatCircleDistance(point1, point2):
	R = 6371e3 #The radius of the Earth in meters
	lat1 = point1[0]
	lat2 = point2[0]
	lon1 = point1[1]
	lon2 = point2[1]
	phi1 = radians(lat1)
	phi2 = radians(lat2)
	deltaPhi = radians(lat2 - lat1)
	deltaLambda = radians(lon2 - lon1)
	a = sin(deltaPhi / 2.0)**2 + cos(phi1)*cos(phi2) * sin(deltaLambda / 2.0)**2
	c = 2 * atan2(sqrt(a), sqrt(1-a))
	d = R * c
	return d

#Reads a file of points where each coordinate is separated by spaces
#and each point is on a new line into a matrix, starting from the column/coordinate
#provided (counting from 0). 
def readFileIntoMatrix(filename, startColumn, rows, cols):
	dataMatrix = np.zeros((rows, cols))
	with open (filename, "r") as dataFile:
		dataString = dataFile.readlines()
		i = 0
		for line in dataString:
			points = line.split()
			j = 0
			for point in points:
				if (j < startColumn):
					pass #do nothing
				else:
					dataMatrix[i][j-startColumn] = point
				j += 1
			i += 1
	return dataMatrix

#Returns the Euclidean distance between two points.
#Points are read as arrays of length numCoords
def euclideanDistance(point1, point2):
	total = 0
	numCoords = len(point1)
	for i in range (0, numCoords):
		total += (point1[i] - point2[i])**2
	return total**(0.5)
		

#Returns the shortest Euclidean distance between any two points in two sets.
#set1 and set2 are arrays of indexes to the points they contain in the dataMatrix
def singleLink(set1, set2, dataMatrix):
	shortestLink = np.inf
	for point in set1:
		for otherPoint in set2:
			currentLink = greatCircleDistance(dataMatrix[point], dataMatrix[otherPoint])
			if (currentLink < shortestLink):
				shortestLink = currentLink
	return shortestLink
	
#Returns the longest Euclidean distance between any two points in two sets.
#set1 and set2 are arrays of indexes to the points they contain in the dataMatrix
def completeLink(set1, set2, dataMatrix):
	longestLink = 0
	for point in set1:
		for otherPoint in set2:
			currentLink = greatCircleDistance(dataMatrix[point], dataMatrix[otherPoint])
			if (currentLink > longestLink):
				longestLink = currentLink
	return longestLink

#Returns a point that represents the average of points in a set
#Points are assumed to have 2 coordinates
def findAverageOfSet(S, dataMatrix):
	xTotal = 0
	yTotal = 0
	for point in S:
		xTotal += dataMatrix[point][0]
		yTotal += dataMatrix[point][1]
	xMean = xTotal / (float)(len(S))
	yMean = yTotal / (float)(len(S))
	return np.array([[xMean], [yMean]])


#Returns the Euclidean distance between the mean position of all points in two sets.
#set1 and set2 are arrays of indexes to the points they contain in the dataMatrix
#assumes points have 2 coordinates (numCoords must be 2)
def meanLink(set1, set2, dataMatrix):
	set1Mean = findAverageOfSet(set1, dataMatrix)
	set2Mean = findAverageOfSet(set2, dataMatrix)
	return greatCircleDistance(set1Mean, set2Mean)
	
def findTwoClosestClusters(sets, dataMatrix, distanceMetric):
	numCoords = dataMatrix.shape[1]
	minDistance = np.inf
	for cluster1 in sets:
		for cluster2 in sets:
			if (cluster1 == cluster2):
				pass
			else:
				distance = distanceMetric(cluster1, cluster2, dataMatrix)
				if (distance < minDistance):
					minDistance = distance
					minClusterPair = [cluster1, cluster2]
	return minClusterPair
	
def runHierarchicalClustering(dataMatrix, distanceMetric, startNumClusters, endNumClusters):
	numCoords = dataMatrix.shape[1]
	numClusters = startNumClusters
	sets = [set() for i in xrange(len(dataMatrix))]
	for i in range (0,len(dataMatrix)):
		sets[i].add(i)
	s1, s2 = findTwoClosestClusters(sets, dataMatrix, distanceMetric)
	while (numClusters > endNumClusters):
		s3 = s1.union(s2)
		sets.remove(s2)
		sets.remove(s1)
		sets.append(s3)
		numClusters -= 1
		s1, s2 = findTwoClosestClusters(sets, dataMatrix, distanceMetric)
	return sets
	
def printClusters(listOfClusteredIndices, dataMatrix):
	k = 1
	for i in listOfClusteredIndices:
		print "Cluster ", k, "is "
		k += 1
		for j in i:
			print dataMatrix[j]
			
#listOfClusteredIndices is a list of sets of indices. Each set corresponds to a cluster, and each index corresponds to a point
#in dataMatrix.
def plotClusters(listOfClusteredIndices, dataMatrix, dotRadius):
	listOfClusters = list()
	k = 0
	for i in listOfClusteredIndices:
		listOfClusters.append(list())
		for j in i:
			listOfClusters[k].append(dataMatrix[j])
		k += 1
	#listOfClusters is a list containing lists of arrays
	for cluster in listOfClusters:
		xCoords = []
		yCoords = []
		for point in cluster:
			xCoords.append(point[0])
			yCoords.append(point[1])
		area = np.pi * (dotRadius)**2
		clusterColor = tuple(np.random.rand(3)) #4 random numbers between 0 and 1
		plt.scatter(xCoords, yCoords, s=area, color=clusterColor)
	plt.show()
	
def gonzalez(dataMatrix, k, initialCenter):
	numPoints = dataMatrix.shape[0]
	setOfCenters = set()
	setOfCenters.add(initialCenter)
	#The indexes of phi are the indexes of points in dataMatrix. The values of phi are the indexes of 
	#the closest center the corresponding point is currently assigned to. 
	phi = np.full(numPoints, initialCenter)
	for i in range(1,k):
		MaxDist = 0
		nextCenter = 0
		for j in range(0, numPoints):
			distToCurrentCenter = greatCircleDistance(dataMatrix[j], dataMatrix[phi[j]])
			if (distToCurrentCenter > MaxDist):
				MaxDist = distToCurrentCenter
				nextCenter = j
		for j in range(0, numPoints):
			distToCurrentCenter = greatCircleDistance(dataMatrix[j], dataMatrix[phi[j]])
			distToNextCenter = greatCircleDistance(dataMatrix[j], dataMatrix[nextCenter])
			if (distToCurrentCenter > distToNextCenter):
				phi[j] = nextCenter
		setOfCenters.add(nextCenter)
	return setOfCenters, phi
	
#Returns a list with a set of point indexes in dataMatrix for each cluster in clusterDict but
#without the key corresponding to the center of the cluster. 
#(Use this for getting sets of points for plotting).
def getListOfClusteredIndicesFromDict(clusterDict):
	listOfClusteredIndices = []
	for k, v in clusterDict.items():
		listOfClusteredIndices.append(v)
	return listOfClusteredIndices
	
#Extracts the cluster information stored in phi from the Gonzalez algorithm.
def getListOfClusteredIndicesFromPhi(phi):
	dictOfClusters = dict()
	for i in range(0, len(phi)):
		key = (int)(phi[i])
		if (key in dictOfClusters):
			dictOfClusters[key].add(i)
		else:
			dictOfClusters[key] = set([i])
	return dictOfClusters, getListOfClusteredIndicesFromDict(dictOfClusters)

	
#Given a list of indices, prints the corresponding points in the data matrix.
def printPointsFromIndices(listOfIndices, dataMatrix):
	for i in listOfIndices:
		print dataMatrix[i]

#Returns the maximum distance between a point in dataMatrix
#and its corresponding center in dictOfClusters
#The keys for dictOfClusters are indices, not actual points		
def getKCenterCost(dictOfClusters, dataMatrix):
	cost = 0
	for k, v in dictOfClusters.items():
		for index in v:
			dist = greatCircleDistance(dataMatrix[index], dataMatrix[k])
			if (dist > cost):
				cost = dist
	return cost
	
#The keys for dictOfClusters are indices, not actual points		
def getKMeansCost(dictOfClusters, dataMatrix):
	total = 0
	for k, v in dictOfClusters.items():
		for index in v:
			total += (greatCircleDistance(dataMatrix[index], dataMatrix[k]))**2
	cost = (total / (float)(dataMatrix.shape[0]))**(0.5)
	return cost
	
#Note: here, the keys of the dictOfClusters are actual points, not indices
def getKMediansCostForPointKeys(dictOfClusters, dataMatrix):
	total = 0
	for k, v in dictOfClusters.items():
		centerPoint = np.array(k)
		for index in v:
			total += (greatCircleDistance(dataMatrix[index], centerPoint))
	cost = total / float(dataMatrix.shape[0])
	return cost
	
#Note: here, the keys of the dictOfClusters are actual points, not indices
def getKMeansCostForPointKeys(dictOfClusters, dataMatrix):
	numCoords = dataMatrix.shape[1]
	total = 0
	for k, v in dictOfClusters.items():
		centerPoint = np.array(k)
		for index in v:
			total += (greatCircleDistance(dataMatrix[index], centerPoint))**2
	cost = (total / (float)(dataMatrix.shape[0]))**(0.5)
	return cost
	
def chooseCWithProbProportionalToDistSquared(dataMatrix, phi):
	#create an array with all the squared distances from x to the c it is closest to.
	squaredDistances = np.zeros([dataMatrix.shape[0], dataMatrix.shape[1]])
	#Find the sum of these distances
	totalSquaredDistance = 0
	for i in range(0, len(phi)):
		squaredDistance = (greatCircleDistance(dataMatrix[i], dataMatrix[phi[i]]))**2
		totalSquaredDistance += squaredDistance
		squaredDistances[i][0] = squaredDistance
		squaredDistances[i][1] = i #to preserve the index number the squared distance in [0] refers to in dataMatrix
		#when we sort the squared distances in ascending order:
	squaredDistances = squaredDistances[squaredDistances[:,1].argsort()]
	#Divide each distance by the sum of distances
	for i in range(0, len(squaredDistances)):
		squaredDistances[i][0] = squaredDistances[i][0]/float(totalSquaredDistance)
	#Make an array with the cumulative sum in each entry
	for i in range(0, len(squaredDistances) - 1):
		squaredDistances[i+1][0] += squaredDistances[i][0]
	#Use a utility to pick a random number, u, between 0 and 1
	u = random.random()
	#Find the first entry in the cumulative density array that u is less than
	#(ie walk from left to right until you find an element that is larger than u)
	#Optionally can put cumulative probability values into binary search tree for faster lookup
	for i in range(0, len(squaredDistances)):
		if (u < squaredDistances[i][0]):
			return squaredDistances[i][1]
	return -1 #error value

	
def kPlusPlus(dataMatrix, k, initialCenter):
	numPoints = dataMatrix.shape[0]
	setOfCenters = set()
	setOfCenters.add(initialCenter)
	#The indexes of phi are the indexes of points in dataMatrix. The values of phi are the indexes of 
	#the closest center the corresponding point is currently assigned to. 
	phi = np.full(numPoints, initialCenter)
	for i in range(1,k):
		#Choose ci from X with probability proportional to d(x, phi_C-1(x))^2
		nextCenter = chooseCWithProbProportionalToDistSquared(dataMatrix, phi)
		#Now update the values in phi for the points that now map to the new center
		for j in range(0, numPoints):
			distToCurrentCenter = greatCircleDistance(dataMatrix[j], dataMatrix[phi[j]])
			distToNextCenter = greatCircleDistance(dataMatrix[j], dataMatrix[nextCenter])
			if (distToCurrentCenter > distToNextCenter):
				phi[j] = nextCenter
		setOfCenters.add(nextCenter)
	return setOfCenters, phi
	

def roundX(x, y):
	return round(x*y)/y

def plotCDF_forKPP(dataMatrix, numTrials, k, startingCenter, roundingFactor):
	threeMeansResults = np.zeros(numTrials)
	maxCost = 0
	totalCost = 0
	minCost = np.inf
	for i in range(0, numTrials):
		kPlusPlusCenters, phi = kPlusPlus(dataMatrix, k, startingCenter)
		dictOfClusters, listOfClusteredIndices = getListOfClusteredIndicesFromPhi(phi)
		cost = getKMeansCost(dictOfClusters, dataMatrix)
		totalCost += cost
		threeMeansResults[i] = cost
		if (cost > maxCost):
			maxCost = cost
		if (cost < minCost):
			minCost = cost
	minCost = roundX(minCost, roundingFactor)
	maxCost = roundX(maxCost, roundingFactor)
	costRange = maxCost - minCost
	bucketSize = 1/float(roundingFactor)
	numBuckets = int(costRange / bucketSize)
	cdf = np.zeros([numBuckets,2])
	j = minCost
	for i in range(0, numBuckets):
		cdf[i][0] = j
		j += bucketSize
	for i in range(0, numTrials):
		roundedCost = roundX(threeMeansResults[i], roundingFactor)
		cdf[roundedCost][1] += 1
	for i in range(0, numBuckets-1):
		cdf[i+1][1] += cdf[i][1]
	for i in range(0, numBuckets):
		cdf[i][1] = cdf[i][1]/float(numTrials)

	fig, ax = plt.subplots(figsize=(8,6))
	plt.plot(cdf[:,0], cdf[:,1])
	ax.grid(True)
	ax.legend(loc='right')
	ax.set_title('CDF of 3-means Cost for K++ Clustering')
	ax.set_xlabel('3-means Cost')
	ax.set_ylabel('Cumulative Density')
	plt.show()
		
#Returns the fraction of the time KPP sorts the data
#into the same subsets (clusters) as Gonzalez
def compareKPPtoGonzalez(dataMatrix, numTrials, k, startingCenter):
	numSame = 0
	for i in range(0, numTrials):
		startingCenter = np.random.randint(0, 1004) #random starting center
		kPlusPlusCenters, phiK = kPlusPlus(dataMatrix, k, startingCenter)
		dictOfClustersK, listOfClusteredIndicesK = getListOfClusteredIndicesFromPhi(phiK)
		gonzalezCenters, phiG = gonzalez(dataMatrix, k, startingCenter)
		dictOfClustersG, listOfClusteredIndicesG = getListOfClusteredIndicesFromPhi(phiG)
		areSame = True
		for cluster in listOfClusteredIndicesK:
			if not cluster in listOfClusteredIndicesG:
				areSame = False
		if areSame:
			numSame +=1
	return numSame / float(numTrials)			

#Compares the difference between keys in two dictionaries
#where the keys are ordered pairs of floats (tuples)
#If keys round to the same number, they are considered the same
def diffBetweenDictKeysIsLarge(dict1, dict2, threshold):
	numCoords = len(dict1.keys()[0])
	setOfKeys1 = set()
	for keyTuple in dict1.keys():
		roundedKey = np.zeros(numCoords)
		for i in range(0, numCoords):
			roundedKey[i] = roundX(keyTuple[i], threshold)
		setOfKeys1.add(tuple(roundedKey))
		
	setOfKeys2 = set()
	for keyTuple in dict2.keys():
		roundedKey = np.zeros(numCoords)
		for i in range(0, numCoords):
			roundedKey[i] = roundX(keyTuple[i], threshold)
		setOfKeys2.add(tuple(roundedKey))
	
	return (setOfKeys1 != setOfKeys2)
	
#map each point to nearest center
def mapEachPointToNearestCenter(dataMatrix, centerTupleDict):
	numPoints = dataMatrix.shape[0]
	numCoords = dataMatrix.shape[1]
	for i in range(0, numPoints):
		minDist = np.inf
		closestCenter = centerTupleDict.keys()[0]#just for starters			
		for center in centerTupleDict.keys():
			centerPoint = np.array(center)
			dist = greatCircleDistance(dataMatrix[i], centerPoint)
			if (dist < minDist):
				minDist = dist
				closestCenter = center
		centerTupleDict[closestCenter].add(i)


def plotCDF_forLloydsUsingKPP(dataMatrix, numTrials, k, roundingFactor, threshold):
	threeMeansResults = np.zeros(numTrials)
	maxCost = 0
	minCost = np.inf
	for i in range(0, numTrials):
		startingCenter = random.randint(0,dataMatrix.shape[0])
		kPlusPlusCenters, phi = kPlusPlus(dataMatrix, k, startingCenter)
		dictOfClusters = lloyd(kPlusPlusCenters, dataMatrix, threshold)
		cost = getKMeansCostForPointKeys(dictOfClusters, dataMatrix)
		threeMeansResults[i] = cost
		if (cost > maxCost):
			maxCost = cost
		if (cost < minCost):
			minCost = cost
	print threeMeansResults
	minCost = roundX(minCost, roundingFactor)
	maxCost = roundX(maxCost, roundingFactor)
	costRange = maxCost - minCost
	if costRange == 0:
		print "All costs were the same, ie:", threeMeansResults[0]
		return
	bucketSize = 1/float(roundingFactor)
	numBuckets = int(costRange / float(bucketSize))
	cdf = np.zeros([numBuckets,2])
	j = minCost
	for i in range(0, numBuckets):
		cdf[i][0] = j
		j += bucketSize
	for i in range(0, numTrials):
		roundedCost = roundX(threeMeansResults[i], roundingFactor)
		cdf[roundedCost][1] += 1
	for i in range(0, numBuckets-1):
		cdf[i+1][1] += cdf[i][1]
	for i in range(0, numBuckets):
		cdf[i][1] = cdf[i][1]/float(numTrials)

	fig, ax = plt.subplots(figsize=(8,6))
	plt.plot(cdf[:,0], cdf[:,1])
	ax.grid(True)
	ax.legend(loc='right')
	ax.set_title('CDF of 3-means Cost for Lloyd\'s Algorithm using start centers from K++')
	ax.set_xlabel('3-means Cost')
	ax.set_ylabel('Cumulative Density')
	plt.show()
		
#Returns the fraction of the time Lloyd's sorts the data
#into the same subsets (clusters) as KPP
def compareLloydToKPP(dataMatrix, numTrials, k, threshold):
	numSame = 0
	for i in range(0, numTrials):
		startingCenter = np.random.randint(0, 1004) #random starting center
		kPlusPlusCenters, phiK = kPlusPlus(dataMatrix, k, startingCenter)
		dictOfClustersK, listOfClusteredIndicesK = getListOfClusteredIndicesFromPhi(phiK)
		dictOfClustersL = lloyd(kPlusPlusCenters, dataMatrix, threshold)
		listOfClusteredIndicesL = getListOfClusteredIndicesFromDict(dictOfClustersL)
		areSame = True
		for cluster in listOfClusteredIndicesK:
			if not cluster in listOfClusteredIndicesL:
				areSame = False
		if areSame:
			numSame +=1
	return numSame / float(numTrials)

#finds the expansion factor needed to increase the radius of a d-dimensional
#ball with radius r so that it has the same volume as the smallest box that
#encloses the ball.
def findExpansionFactorForBallRadius(d):
	return 2 * ((math.gamma(d/float(2) + 1))**(1/float(d))/(math.pi)**(0.5))
	
#Plots the expansion factor found above from d=zero-ish to d=d-1
def plotExpansionFactor(d):
	expFactors = np.zeros(d)
	expFactors[0] = findExpansionFactorForBallRadius(0.001) # to avoid division by zero
	for i in range(1,d):
		expFactors[i] = findExpansionFactorForBallRadius(i)
		
	fig, ax = plt.subplots()
	ax.set_title('High-Dimensional Distances')
	ax.set_xlabel('dimension d')
	ax.set_ylabel('expansion factor c')
	plt.plot(expFactors)
	plt.show()
	
def lloyd(centers, dataMatrix, threshold):
	#First convert centers from indices in dataMatrix into tuples with the
	#actual values in them. Make each tuple the key in a dict, referring
	#to an empty set, where the indices of the points in the cluster will go.
	numPoints = dataMatrix.shape[0]
	numCoords = dataMatrix.shape[1]
	centerTupleDict = dict()
	for i in centers:
		thisTuple = tuple(dataMatrix[i])
		centerTupleDict[thisTuple] = set()
		
	newCenterTupleDict = dict()#garbage initialization value
	while diffBetweenDictKeysIsLarge(centerTupleDict, newCenterTupleDict, threshold):
		mapEachPointToNearestCenter(dataMatrix, centerTupleDict)
		#Now reassign every center to the average of the points that are assigned to it
		newCenterTupleDict.clear() #We only wanted newCenterTupleDict's keys for the while loop comparison;
		#its sets are outdated now.
		for center, pointSet in centerTupleDict.items():
			totals = np.zeros(numCoords)
			averages = np.zeros(numCoords)
			for pointIndex in pointSet:
				for i in range(0, numCoords):
					totals[i] += dataMatrix[pointIndex][i]
				for i in range(0, numCoords):
					averages[i] = totals[i] / float(len(pointSet))
			newTuple = tuple(averages)
			newCenterTupleDict[newTuple] = set()
			
		tempDict = newCenterTupleDict
		newCenterTupleDict = centerTupleDict
		centerTupleDict = tempDict
	
	#Do this mapping one last time to get most up-to-date clusters
	mapEachPointToNearestCenter(dataMatrix, centerTupleDict)	
	return centerTupleDict
	
def getPointsInOneDimensionFromSetOfIndices(setOfIndices, dataMatrix, dimensionIndex):
	numPoints = len(setOfIndices)
	points = np.zeros(numPoints)
	i = 0
	for index in setOfIndices:
		points[i] = dataMatrix[index][dimensionIndex]
		i += 1
	return points
	
	
def kMediansClustering(centers, dataMatrix, threshold):
	#First convert centers from indices in dataMatrix into tuples with the
	#actual values in them. Make each tuple the key in a dict, referring
	#to an empty set, where the indices of the points in the cluster will go.
	k = len(centers)
	numPoints = dataMatrix.shape[0]
	numCoords = dataMatrix.shape[1]
	centerTupleDict = dict()
	for i in centers:
		thisTuple = tuple(dataMatrix[i])
		centerTupleDict[thisTuple] = set()
		
	newCenterTupleDict = dict()#garbage initialization value
	while diffBetweenDictKeysIsLarge(centerTupleDict, newCenterTupleDict, threshold):
		mapEachPointToNearestCenter(dataMatrix, centerTupleDict)
		#Now reassign every center to the median in each dimension of the points that are assigned to it
		newCenterTupleDict.clear() #We only wanted newCenterTupleDict's keys for the while loop comparison;
		#its sets are outdated now.
		for center, pointSet in centerTupleDict.items():
			median = np.zeros(numCoords)
			for i in range(0, numCoords):
				points = getPointsInOneDimensionFromSetOfIndices(pointSet, dataMatrix, i)
				sortedPoints = np.sort(points)
				median[i] = sortedPoints[len(sortedPoints) / 2]
			newCenterTupleDict[tuple(median)] = set()
		tempDict = newCenterTupleDict
		newCenterTupleDict = centerTupleDict
		centerTupleDict = tempDict
	
	#Do this mapping one last time to get most up-to-date clusters
	mapEachPointToNearestCenter(dataMatrix, centerTupleDict)	
	return centerTupleDict


def combineMatrices(AllMatricesList):
	totalRows = 0
	for matrix in AllMatricesList:
		totalRows += len(matrix)
	combinedMatrix = np.zeros([totalRows, 2])
	for matrix in AllMatricesList:
		i = 0
		for row in matrix:
			combinedMatrix[i][0] = row[0]
			combinedMatrix[i][1] = row[1]
			i += 1
	return combinedMatrix

def plotCostVersusK(dataMatrix, clusteringAlgorithm, costFunction, minK, maxK, stepSize):
	costArray = np.zeros(maxK)
	initialCenter = np.random.randint(0,len(dataMatrix))
	for i in range(minK, maxK, stepSize):
		centers, phi = clusteringAlgorithm(dataMatrix, i, initialCenter)
		dictOfClusters, listOfClusteredIndices = getListOfClusteredIndicesFromPhi(phi)
		cost = costFunction(dictOfClusters, dataMatrix)
		costArray[i] = cost
		if i == minK:
			pass
		else: #We fill in the intervening values in the array with a linear function
		#to approximate the intervening cost values
			previousCost = costArray[i - stepSize]
			slope = (cost - previousCost) / (stepSize)
			intercept = cost - slope * i
			for j in range((i - stepSize)+1, i):
				costArray[j] = intercept + j * slope
	fig, ax = plt.subplots()
	ax.set_title('Cost versus k')
	ax.set_xlabel('k')
	ax.set_ylabel('Cost')
	plt.plot(costArray)
	plt.show()



		

#"main" section:
ReptilesMatrix = readCSVFileIntoMatrix("Reptiles.txt", 27, 28)
TerrestrialMammalsMatrix = readCSVFileIntoMatrix("TerrestrialMammals.txt", 27, 28)
AmphibiansMatrix = readCSVFileIntoMatrix("Amphibians.txt", 27, 28)
AllMatricesList = list([ReptilesMatrix, TerrestrialMammalsMatrix, AmphibiansMatrix])
AllAnimalsMatrix = combineMatrices(AllMatricesList)


#plotCostVersusK(AllAnimalsMatrix, gonzalez, getKMeansCost, 0, 11, 5)

#point1 = np.array([0,180])
#point2 = np.array([0,0])
#print greatCircleDistance(point1, point2) / 1000.0


#singleLinkClusters = runHierarchicalClustering(ReptilesMatrix, singleLink, ReptilesMatrix.shape[0], 10)
###printClusters(singleLinkClusters, C1DataMatrix)
#plotClusters(singleLinkClusters, ReptilesMatrix, 10)

k=40
dotRadius = 0.3
initialCenter = np.random.randint(0,len(AllAnimalsMatrix))
#gonzalezCenters, phi = gonzalez(AllAnimalsMatrix, k, initialCenter)
###printPointsFromIndices(gonzalezCenters, AllAnimalsMatrix)
#dictOfClusters, listOfClusteredIndices = getListOfClusteredIndicesFromPhi(phi)
#plotClusters(listOfClusteredIndices, AllAnimalsMatrix, dotRadius)
###print getKCenterCost(dictOfClusters, )
#print "The centers are", printPointsFromIndices(gonzalezCenters, AllAnimalsMatrix)
#print "Kmeans cost is", getKMeansCost(dictOfClusters, AllAnimalsMatrix)
#


#kPlusPlusCenters, phi = kPlusPlus(AllAnimalsMatrix, k, initialCenter)
##printPointsFromIndices(kPlusPlusCenters, AllAnimalsMatrix)
#dictOfClusters, listOfClusteredIndices = getListOfClusteredIndicesFromPhi(phi)
#plotClusters(listOfClusteredIndices, AllAnimalsMatrix, dotRadius)

#print getKCenterCost(dictOfClusters, ReptilesMatrix)
#print getKMeansCost(dictOfClusters, ReptilesMatrix)
#numTrials = 20
#k = 3
#roundingFactor = 128
#plotCDF_forKPP(ReptilesMatrix, numTrials, k, initialCenter, roundingFactor)
#
#print compareKPPtoGonzalez(ReptilesMatrix, numTrials, k, initialCenter)

k=40
threshold = 128
#startingCenters = list()
#for i in range (0,k):
#	startingCenters.append(i)
#startingCenters = set(startingCenters)
#print startingCenters
startingCenters, phi = gonzalez(AllAnimalsMatrix, k, initialCenter)
dictOfClusters = lloyd(startingCenters, AllAnimalsMatrix, threshold)
listOfClusters = getListOfClusteredIndicesFromDict(dictOfClusters)
plotClusters(listOfClusters, AllAnimalsMatrix, dotRadius)
#print getKMeansCostForPointKeys(dictOfClusters, ReptilesMatrix)
#print dictOfClusters.keys()









#dotRadius = 2
#threshold = 128
#startingCenters = set([0,1,2])
#dictOfClusters = lloyd(startingCenters, ReptilesMatrix, threshold)
#listOfClusters = getListOfClusteredIndicesFromDict(dictOfClusters)
#plotClusters(listOfClusters, ReptilesMatrix, dotRadius)
#print getKMeansCostForPointKeys(dictOfClusters, ReptilesMatrix)
#print dictOfClusters.keys()
#
#k=3
#initialCenter = 0
#threshold = 128
#dotRadius = 2
#numTrials = 20
#roundingFactor = 128
#gonzalezCenters, phi = gonzalez(ReptilesMatrix, k, initialCenter)
#dictOfClusters = lloyd(gonzalezCenters, ReptilesMatrix, threshold)
#listOfClusters = getListOfClusteredIndicesFromDict(dictOfClusters)
#plotClusters(listOfClusters, ReptilesMatrix, dotRadius)
#print getKMeansCostForPointKeys(dictOfClusters, ReptilesMatrix)
#print dictOfClusters.keys()
#plotCDF_forLloydsUsingKPP(ReptilesMatrix, numTrials, k, roundingFactor, threshold)
#print compareLloydToKPP(ReptilesMatrix, numTrials, k, threshold)
#
#print findExpansionFactorForBallRadius(2)
#print findExpansionFactorForBallRadius(3)
#print findExpansionFactorForBallRadius(4)
#plotExpansionFactor(21)
#
#k=4
#threshold = 128
#startingCenters = set([0,1,2,3])
#dictOfClustersLloyd = lloyd(startingCenters, C3DataMatrix, threshold)
#dictOfClustersKMedians = kMediansClustering(startingCenters, C3DataMatrix, threshold)
#print "4 Means cost for Lloyd:",getKMeansCostForPointKeys(dictOfClustersLloyd, C3DataMatrix)
#print "4 Means cost for kMedians:",getKMeansCostForPointKeys(dictOfClustersKMedians, C3DataMatrix)
#print "4 Medians cost for Lloyd:",getKMediansCostForPointKeys(dictOfClustersLloyd, C3DataMatrix)
#print "4 Medians cost for kMedians:",getKMediansCostForPointKeys(dictOfClustersKMedians, C3DataMatrix)
#
#gonzalezCenters, phi = gonzalez(C3DataMatrix, k, initialCenter)
#dictOfClustersLloyd = lloyd(gonzalezCenters, C3DataMatrix, threshold)
#dictOfClustersKMedians = kMediansClustering(gonzalezCenters, C3DataMatrix, threshold)
#print "4 Means cost for Lloyd:",getKMeansCostForPointKeys(dictOfClustersLloyd, C3DataMatrix)
#print "4 Means cost for kMedians:",getKMeansCostForPointKeys(dictOfClustersKMedians, C3DataMatrix)
#print "4 Medians cost for Lloyd:",getKMediansCostForPointKeys(dictOfClustersLloyd, C3DataMatrix)
#print "4 Medians cost for kMedians:",getKMediansCostForPointKeys(dictOfClustersKMedians, C3DataMatrix)
#print "kMedians centers:", dictOfClustersKMedians.keys()
#
#k=5
#threshold = 128
#initialCenter=0
#gonzalezCenters, phi = gonzalez(C3DataMatrix, k, initialCenter)
#dictOfClustersLloyd = lloyd(gonzalezCenters, C3DataMatrix, threshold)
#dictOfClustersKMedians = kMediansClustering(gonzalezCenters, C3DataMatrix, threshold)
#print "5 Means cost for Lloyd:",getKMeansCostForPointKeys(dictOfClustersLloyd, C3DataMatrix)
#print "5 Means cost for kMedians:",getKMeansCostForPointKeys(dictOfClustersKMedians, C3DataMatrix)
#print "5 Medians cost for Lloyd:",getKMediansCostForPointKeys(dictOfClustersLloyd, C3DataMatrix)
#print "5 Medians cost for kMedians:",getKMediansCostForPointKeys(dictOfClustersKMedians, C3DataMatrix)
#print "kMedians centers:", dictOfClustersKMedians.keys()

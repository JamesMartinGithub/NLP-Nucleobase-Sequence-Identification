import math
import numpy as np
import random
from sklearn.cluster import KMeans
from scipy import spatial

filefolders = ["Human", "Rat"]
filenames = ["ATR", "BUB3", "CDC20", "CDK2", "ERH", "NEK2", "NUMA1", "PIN1", "PRC1", "SGO1"]
seenTestFileNames = ["SGO1_H", "SGO1_R", "NUMA1_H", "NUMA1_R", "NEK2_H", "NEK2_R", "CDC20_H", "CDC20_R", "CDK2_H", "CDK2_R", "ERH_H", "ERH_R", "ATR_H", "ATR_R", "BUB3_H", "BUB3_R", "PIN1_H", "PIN1_R", "PRC1_H", "PRC1_R"]
unseenTestFileNames = ["CDK1_H", "CDK1_R", "CDKN1B_H", "CDKN1B_R", "FST_H", "FST_R", "KHDRBS1_H", "KHDRBS1_R", "NEK4_H", "NEK4_R", "NEK8_H", "NEK8_R", "PLK1_H", "PLK1_R", "RALBP1_H", "RALBP1_R"]
clusterNum = 200  # Number of k-means clusters
sampleNum = 8000  # Number of samples used in k-means model fitting, randomly selected
vecLength = 10  # Length of vectors taken from nucleotide sequence
nRepeat = 2  # Number of repeated nucleotides that forms a vector

histograms = []
kmeans = KMeans(clusterNum)

def train():
    global histograms
    # Get all vectors
    genes = []
    joinedVecs = []
    for folder in filefolders:
        for name in filenames:
            file = open("Train\\" + folder + "\\" + name + ".fna", "r")
            file.readline()
            geneVecs = []
            for line in file:
                geneVecs.extend(fileToVecs(line))
            # Add gene's vector list to gene list
            genes.append(geneVecs)
            joinedVecs.extend(geneVecs)
    # K-means to quantise
    kmeans.fit_predict(select(sampleNum, joinedVecs))
    # Get labels
    labels = kmeans.predict(joinedVecs)
    # Use labels to create histograms for each gene
    histograms = createHistograms(genes, labels)
    print("Trained")


def fileToVecs(line):
    geneVecs = []
    if line[0] == '>' or len(line[:-1]) < (nRepeat + vecLength):
        return geneVecs
    else:
        # Create vectors and add to gene's vector list
        seq = line[:-1]
        sLength = nRepeat + vecLength
        for i in range(0, len(seq) - (nRepeat + vecLength)):
            # Vectors are taken after double repeated nucleotides for shift invariance
            if seq[i:i+nRepeat] == "A" * nRepeat and seq[i+nRepeat] != "A":
                geneVecs.append(seqToVec(seq[i+nRepeat:i+sLength]))
            if seq[i:i+nRepeat] == "T" * nRepeat and seq[i+nRepeat] != "T":
                geneVecs.append(seqToVec(seq[i+nRepeat:i+sLength]))
            if seq[i:i+nRepeat] == "C" * nRepeat and seq[i+nRepeat] != "C":
                geneVecs.append(seqToVec(seq[i+nRepeat:i+sLength]))
            if seq[i:i+nRepeat] == "G" * nRepeat and seq[i+nRepeat] != "G":
                geneVecs.append(seqToVec(seq[i+nRepeat:i+sLength]))
    return geneVecs


# Convert nucleotide sequence string to vector of length 10
def seqToVec(seq):
    vector = np.zeros(vecLength, dtype=float)
    for i in range(0, vecLength):
        if seq[i] == 'A':
            vector[i] = 1
        if seq[i] == 'T':
            vector[i] = 5
        if seq[i] == 'C':
            vector[i] = 10
        if seq[i] == 'G':
            vector[i] = 15
    size = 0
    for v in vector:
        size += pow(v, 2)
    vector /= math.sqrt(size)
    return vector


# Randomly selects vectors from list of vectors, given number to select
def select(number, vectors):
    selected = []
    indexes = list(range(len(vectors)))
    random.shuffle(indexes)
    for i in range(number):
        selected.append(vectors[indexes[i]])
    return np.array(selected)


# Returns list of histograms of word counts given list of quantised vectors, per category
def createHistograms(categoryVectors, labels):
    histogramList = np.zeros((20, clusterNum))
    startIndex = 0
    for x in range(0, 20):
        for y in range(startIndex, startIndex + len(categoryVectors[x])):
            if labels[y] != -1:
                histogramList[x][labels[y]] += 1
        startIndex += len(categoryVectors[x])
    return histogramList


# Returns prediction accuracy, calculated by comparing predicted category to actual category
def test():
    def testOne(name, seen: bool):
        # Get test vectors
        if seen:
            file = open("Test\\Seen\\" + name + ".fa", "r")
        else:
            file = open("Test\\Unseen\\" + name + ".fna", "r")
        file.readline()
        geneVecs = []
        for line in file:
            geneVecs.extend(fileToVecs(line))
        # Get labels
        labels = kmeans.predict(geneVecs)
        # Create histogram
        histogram = createHistogram(labels)
        # Classify histogram
        category, gene = classify(histogram, histograms)
        # Return result
        #print(name + " is " + category)
        categoryInt, geneInt = (0, 0)
        if category[0] == name[len(name) - 1]:
            categoryInt = 1
        else:
            categoryInt = 0
        if gene == name[:-2] and categoryInt == 1:
            geneInt = 1
        else:
            geneInt = 0
        return (categoryInt, geneInt)
    seenCorrect, seenGeneCorrect, unseenCorrect = (0, 0, 0)
    # Test each test file
    for name in seenTestFileNames:
        (cat, gen) = testOne(name, True)
        seenCorrect += cat
        seenGeneCorrect += gen
    for name in unseenTestFileNames:
        unseenCorrect += testOne(name, False)[0]
    # Calculate accuracies
    seenAccuracy = ((seenCorrect / len(seenTestFileNames)) * 100)
    unseenAccuracy = ((unseenCorrect / len(unseenTestFileNames)) * 100)
    seenGeneAccuracy = ((seenGeneCorrect / len(seenTestFileNames)) * 100)
    #print("Seen Accuracy:", str(seenAccuracy) + "%")
    #print("Unseen Accuracy:", str(unseenAccuracy) + "%")
    return (seenAccuracy, unseenAccuracy, seenGeneAccuracy)


# Returns histogram of word counts given list of quantised vectors
def createHistogram(labels):
    histogram = np.zeros(clusterNum)
    for i in range(len(labels)):
        if labels[i] != -1:
            histogram[labels[i]] += 1
    return histogram


# Finds closest matching gene histogram to snippet histogram, to return closest category
def classify(histogram, histogramList):
    highScore = 0
    highestCategory = -1
    tupleList = []
    # Find the best similarity between snippet histogram and all gene histograms
    for i in range(len(histogramList)):
        # Calculate cosine similarity between histograms
        similarity = 1 - spatial.distance.cosine(histogram, histogramList[i])
        if i <= 9:
            tupleList.append((filenames[i % 10] + "_H", similarity))
        else:
            tupleList.append((filenames[i % 10] + "_R", similarity))
        if similarity > highScore:
            highScore = similarity
            highestCategory = i
    #print("Gene: " + filenames[highestCategory % 10])
    #print(sorted(tupleList, key=lambda x: x[1], reverse=True))
    if highestCategory <= 9:
        return (filefolders[0], filenames[highestCategory % 10])
    else:
        return (filefolders[1], filenames[highestCategory % 10])

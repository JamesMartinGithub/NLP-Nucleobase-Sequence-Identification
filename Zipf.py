import math
import numpy as np
from scipy import spatial

filefolders = ["Human", "Rat"]
filenames = ["ATR", "BUB3", "CDC20", "CDK2", "ERH", "NEK2", "NUMA1", "PIN1", "PRC1", "SGO1"]
seenTestFileNames = ["SGO1_H", "SGO1_R", "NUMA1_H", "NUMA1_R", "NEK2_H", "NEK2_R", "CDC20_H", "CDC20_R", "CDK2_H", "CDK2_R", "ERH_H", "ERH_R", "ATR_H", "ATR_R", "BUB3_H", "BUB3_R", "PIN1_H", "PIN1_R", "PRC1_H", "PRC1_R"]
unseenTestFileNames = ["CDK1_H", "CDK1_R", "CDKN1B_H", "CDKN1B_R", "FST_H", "FST_R", "KHDRBS1_H", "KHDRBS1_R", "NEK4_H", "NEK4_R", "NEK8_H", "NEK8_R", "PLK1_H", "PLK1_R", "RALBP1_H", "RALBP1_R"]
permutations = ["AAAA", "AATT", "AACC", "AAGG", "TTAA", "TTTT", "TTCC", "TTGG", "CCAA", "CCTT", "CCCC", "CCGG", "GGAA", "GGTT", "GGCC", "GGGG"]

histograms = np.zeros((20, len(permutations)))


def train():
    global histograms
    # Import training nucleotide sequences into list of histograms of base counts
    index = 0
    for folder in filefolders:
        for name in filenames:
            file = open("Train\\" + folder + "\\" + name + ".fna", "r")
            file.readline()
            for line in file:
                for i in range(0, len(line) - (1 + len(permutations[0]))):
                    i, value = histChangeValue(line[i:i+len(permutations[0])])
                    histograms[index, i] += value
            # Normalise histogram
            size = 0
            for v in histograms[index]:
                size += pow(v, 2)
            histograms[index] /= math.sqrt(size)
            index += 1


# Return index of variable to change, and value to change variable by
def histChangeValue(c):
    if c in permutations:
        return permutations.index(c), 1
    else:
        return 0, 0

def test():
    def testOne(name, seen: bool):
        # Get test vectors
        if seen:
            file = open("Test\\Seen\\" + name + ".fa", "r")
        else:
            file = open("Test\\Unseen\\" + name + ".fna", "r")
        file.readline()
        # Create histogram of test data
        histogram = np.zeros(len(permutations))
        for line in file:
            for i in range(0, len(line) - (1 + len(permutations[0]))):
                i, value = histChangeValue(line[i:i + len(permutations[0])])
                histogram[i] += value
        #Normalise histogram
        size = 0
        for v in histogram:
            size += pow(v, 2)
        histogram /= math.sqrt(size)
        # Classify histogram to the closest trained histogram
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
    #train()
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
    return (seenAccuracy, unseenAccuracy, seenGeneAccuracy)


# Finds closest matching gene histogram to test histogram, to return closest category
def classify(histogram, histogramList):
    highScore = 0
    highestCategory = -1
    tupleList = []
    # Find the best similarity between test histogram and all gene histograms
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

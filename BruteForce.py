filefolders = ["Human", "Rat"]
filenames = ["ATR", "BUB3", "CDC20", "CDK2", "ERH", "NEK2", "NUMA1", "PIN1", "PRC1", "SGO1"]
seenTestFileNames = ["SGO1_H", "SGO1_R", "NUMA1_H", "NUMA1_R", "NEK2_H", "NEK2_R", "CDC20_H", "CDC20_R", "CDK2_H", "CDK2_R", "ERH_H", "ERH_R", "ATR_H", "ATR_R", "BUB3_H", "BUB3_R", "PIN1_H", "PIN1_R", "PRC1_H", "PRC1_R"]
unseenTestFileNames = ["CDK1_H", "CDK1_R", "CDKN1B_H", "CDKN1B_R", "FST_H", "FST_R", "KHDRBS1_H", "KHDRBS1_R", "NEK4_H", "NEK4_R", "NEK8_H", "NEK8_R", "PLK1_H", "PLK1_R", "RALBP1_H", "RALBP1_R"]
checkLength = 10
checkScore = 0.6

sequences = [""] * 20


def train():
    global sequences
    # Import training nucleotide sequences into list of strings
    index = 0
    for folder in filefolders:
        for name in filenames:
            file = open("Train\\" + folder + "\\" + name + ".fna", "r")
            file.readline()
            for line in file:
                sequences[index] += line[:-1]
            index += 1

def test():
    def testOne(name, seen: bool):
        # Import test data
        if seen:
            file = open("Test\\Seen\\" + name + ".fa", "r")
        else:
            file = open("Test\\Unseen\\" + name + ".fna", "r")
        file.readline()
        testSequence = ""
        for line in file:
            testSequence += line[:-1]
        # Check each known gene nucleotide sequence for test sequence
        highestScore = 0
        bestCategory = "H"
        index = 0
        bestIndex = 0
        categoryInt, geneInt = (0, 0)
        for seq in sequences:
            # Check every sub-sequence
            for i in range(0, len(seq) - checkLength):
                # Check if sub-sequence is worth checking by looking at first 'checkLength' characters
                if scoreMatch(testSequence, seq[i:i + checkLength], limit=checkLength) > checkScore:
                    score = scoreMatch(testSequence, seq[i:])
                    if score > highestScore:
                        highestScore = score
                        if index <= 9:
                            bestCategory = "H"
                        else:
                            bestCategory = "R"
                        bestIndex = index
            index += 1
        #print(name, "=", filenames[bestIndex % 10])
        if bestCategory == name[len(name) - 1]:
            categoryInt = 1
        else:
            categoryInt = 0
        if filenames[bestIndex % 10] == name[:-2] and categoryInt == 1:
            geneInt = 1
        else:
            geneInt = 0
        return (categoryInt, geneInt)
    train()
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


# Return score 0.0-1.0 of match between test sequence and known sequence
def scoreMatch(testSeq, trainSeq, limit: int = 0):
    score = 0.0
    charCount = 0
    testLength = len(testSeq)
    for c in trainSeq:
        if (limit != 0 and charCount >= limit) or (charCount >= testLength):
            break
        if c == testSeq[charCount]:
            score += 1
        charCount += 1
    return score / charCount

import BagOfWords
import BruteForce
import Zipf
import time

repeats = 20


# Bag of Words
def bow():
    seenSum, unseenSum, seenGeneSum, averageTime = (0, 0, 0, 0.0)
    for i in range(0, repeats):
        BagOfWords.train()
        t0 = time.time()
        seenCorrect, unseenCorrect, seenGeneCorrect = BagOfWords.test()
        t1 = time.time()
        averageTime += t1 - t0
        seenSum += seenCorrect
        unseenSum += unseenCorrect
        seenGeneSum += seenGeneCorrect
    print("Average time:", averageTime / repeats)
    print("Seen average accuracy:", str(seenSum / repeats) + "%")
    print("  Seen gene average accuracy:", str(seenGeneSum / repeats) + "%")
    print("Unseen average accuracy:", str(unseenSum / repeats) + "%")
    print("Total average accuracy:", str((seenSum + unseenSum) / (repeats * 2)) + "%")


#Brute Force
def brute():
    seenSum, unseenSum, seenGeneSum, averageTime = (0, 0, 0, 0)
    for i in range(0, repeats):
        t0 = time.time()
        seenCorrect, unseenCorrect, seenGeneCorrect = BruteForce.test()
        t1 = time.time()
        averageTime += t1 - t0
        seenSum += seenCorrect
        unseenSum += unseenCorrect
        seenGeneSum += seenGeneCorrect
    print("Average time:", averageTime / repeats)
    print("Seen average accuracy:", str(seenSum / repeats) + "%")
    print("  Seen gene average accuracy:", str(seenGeneSum / repeats) + "%")
    print("Unseen average accuracy:", str(unseenSum / repeats) + "%")
    print("Total average accuracy:", str((seenSum + unseenSum) / (repeats * 2)) + "%")


# Frequency Analysis
def zipf():
    seenSum, unseenSum, seenGeneSum, averageTime = (0, 0, 0, 0)
    for i in range(0, repeats):
        Zipf.train()
        t0 = time.time()
        seenCorrect, unseenCorrect, seenGeneCorrect = Zipf.test()
        t1 = time.time()
        averageTime += t1 - t0
        seenSum += seenCorrect
        unseenSum += unseenCorrect
        seenGeneSum += seenGeneCorrect
    print("Average time:", averageTime / repeats)
    print("Seen average accuracy:", str(seenSum / repeats) + "%")
    print("  Seen gene average accuracy:", str(seenGeneSum / repeats) + "%")
    print("Unseen average accuracy:", str(unseenSum / repeats) + "%")
    print("Total average accuracy:", str((seenSum + unseenSum) / (repeats * 2)) + "%")


if __name__ == '__main__':
    # Bag of Words
    bow()

    # Brute Forces
    #brute()

    # Zipf
    #zipf()

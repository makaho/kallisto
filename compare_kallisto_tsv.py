#!/usr/bin/python
import os.path
import csv
import sys


#call with count file
def laodCounts(countFilename):
    counts = dict()
    with open(countFilename) as csv_counts:
        reader_csv_counts = csv.reader(csv_counts, delimiter='\t')
        for row in reader_csv_counts:
            counts[int(row[0])] = int(row[1])
    return counts

#call with ECs file
def parseECfile(ecFilename):
    multiEcInverseIndex = dict()
    with open(ecFilename) as csv_counts_a:
        reader_csv_counts_a = csv.reader(csv_counts_a, delimiter='\t')
        count = 1
        for row in reader_csv_counts_a:
            if row[0] == row[1]:
                count += 1
            else:
                multiEcInverseIndex[row[1]] = row[0]
    num_lines = sum(1 for line in open(ecFilename))
    return num_lines, count, multiEcInverseIndex


def compare(dirA, dirB):
    dirA = os.path.join(dirA, '')
    dirB = os.path.join(dirB, '')
    counts_a = laodCounts(dirA + "pseudoalignments.tsv")
    counts_b = laodCounts(dirB + "pseudoalignments.tsv")

    total_a, oneECs_a, inverseECs_a = parseECfile(dirA + "pseudoalignments.ec")
    total_b, oneECs_b, inverseECs_b = parseECfile(dirB + "pseudoalignments.ec")

    error = False

    if total_a != total_b:
        print("total number of ECs not equal")
        print("a: " + str(total_a))
        print("b: " + str(total_b))
        error = True
    if oneECs_a != oneECs_b:
        print("number of oneECs not equal")
        print("a: " + str(oneECs_a))
        print("b: " + str(oneECs_b))
        error = True
    if set(inverseECs_a.keys()) != set(inverseECs_b.keys()):
        print("multi-EC classes are not equal")
        print(set(inverseECs_a.keys()).symmetric_difference(set(inverseECs_b.keys())))
        print(len(set(inverseECs_a.keys()).symmetric_difference(set(inverseECs_b.keys()))))
        error = True
    for i in range(0, oneECs_a):
        if int(counts_a[i]) != int(counts_b[i]):
            if error == False:
                print("count not equal")
                print("for single-EC key " + str(i))
                print(str(counts_a[i]))
                print(str(counts_b[i]))
                error = True
    for key in inverseECs_a.keys():
        if counts_a[int(inverseECs_a[key])] != counts_b[int(inverseECs_b[key])]:
            if error == False:
                print("count not equal")
                print("for multi-EC key " + key)
                print("EC in a " + str(inverseECs_a[key]))
                print("EC in b " + str(inverseECs_b[key]))
                print("value in a " + str(counts_a[int(inverseECs_a[key])]))
                print("value in b " + str(counts_b[int(inverseECs_b[key])]))
                error = True
    if error == True:
        print("Files are NOT equal :(")
    if error == False:
        print("Files are equal :)")
    return error

def main():
    # count args
    if len(sys.argv) < 3:
        print("usage " + sys.argv[0] + " directory1 directory2")

    dirA = sys.argv[1];
    dirB = sys.argv[2];
    compare(dirA, dirB)


if __name__ == "__main__":
    main()

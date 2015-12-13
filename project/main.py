import fileinput
import collections
from collections import defaultdict
import os
import csv
import string
import random

map = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"STOP", "TAG":"STOP",
    "TGT":"C", "TGC":"C", "TGA":"STOP", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

def parseGenomeFile(f):
    data = collections.defaultdict(list)
    key = ''
    val = ''
    first = True
    for line in f:
        if '>' in line:
            if not first:
                data[key] = val
                val = ''
            key = line.rstrip()
            first = False
        else:
            val += line.rstrip()
    else:
        data[key] = val

    return data

def setupDB():
    buildDB = True
    ''' Find files to make up database '''
    dbFiles = set()
    global db
    files = set()
    fileSize = collections.defaultdict(list)
    for f in os.listdir('DB_Files/'):
        if f.endswith('.fa') or f.endswith('.fsa') or f.endswith('.fasta'):
            dbFiles.add(f.rstrip())
    for f in os.listdir('.'):
        if f == 'db.csv':
            for line in file('db.csv'):
                if 'File,' in line:
                    t = line[line.find(',')+1: line.rfind(',')]
                    files.add(t)
                    fileSize[t] = line[line.rfind(',')+1:].rstrip()

    for item in files:
        if item not in dbFiles:
            buildDB = True
            break
        else:
            for i in fileSize.items():
                valI = i[1].rstrip()
                valItem = os.path.getsize('DB_Files/' + item)
                if i[0] == item:
                    if int(valI) == int(valItem):
                        buildDB = False
                    else:
                        buildDB = True
                        break
                if buildDB == True:
                    break

    ''' Read in db.csv if it exists '''
    if not buildDB:
        first = True
        read = False
        with open('db.csv') as f:
            for i in xrange(len(dbFiles)+1):
                f.next()
            for line in f:
                key = line[:line.find(',')].rstrip()
                val = line[line.find(',')+1:].rstrip()
                db[key] = val
        return db

    ''' Generate db from the set of files'''
    for f in dbFiles:
        f = file('DB_Files/' + f)
        key = ''
        val = ''
        first = True
        for line in f:
            if '>' in line:
                if not first:
                    db[val] = key
                    val = ''
                key = line.rstrip()
                first = False
            else:
                val += line.rstrip()
        else:
            db[val] = key
            val = ''
    else:
        ''' Saves database to db.csv file '''
        w = csv.writer(open("db.csv", "w"))
        w.writerow(['Length', len(dbFiles)])
        for f in dbFiles:
            w.writerow(['File', f, os.path.getsize('DB_Files/' + f)])
        for key, val in db.items():
            w.writerow([key, val])
    return db

''' Taken from Dr. Langmead HW5 solutions '''
_revcomp_trans = string.maketrans("ACGTacgt", "TGCAtgca")
def getComplementStrand(s):
    return s[::-1].translate(_revcomp_trans)

def getProteinData(strand):
    proteins = set()
    startCodon = 'ATG'
    stop = ('TAA', 'TGA', 'TAG')
    global db
    
    ''' Modified from Dr. Langmead's HW5 solution '''
    for seq in [strand, getComplementStrand(strand)]:
        for x in xrange(0, len(seq), 1): 
            codon = seq[x:x+3]
            if len(codon) < 3: break
            if codon == startCodon:
                proteinSeq = ''
                for i in xrange(x, len(seq), 3):
                    codon2 = seq[i:i+3]
                    if len(codon2) < 3: break
                    if codon2 in stop:
                        proteins.add(proteinSeq)
                        break
                    temp = random.random()
                    if codon2 in map.keys():
                        proteinSeq += map[codon2]
                    elif 'K' in codon2:
                        if temp < 0.5:
                            proteinSeq += map[codon2.replace('K', 'G')]
                        else:
                            proteinSeq += map[codon2.replace('K', 'T')]
                    elif 'M' in codon2:
                        if temp < 0.5:
                            proteinSeq += map[codon2.replace('M', 'A')]
                        else:
                            proteinSeq += map[codon2.replace('M', 'C')]
                    elif 'B' in codon2:
                        if temp < 0.33:
                            proteinSeq += map[codon2.replace('B', 'G')]
                        elif temp > 0.67:
                            proteinSeq += map[codon2.replace('B', 'C')]
                        else:
                            proteinSeq += map[codon2.replace('B', 'T')]
                    elif 'V' in codon2:
                        if temp < 0.33:
                            proteinSeq += map[codon2.replace('V', 'G')]
                        elif temp > 0.67:
                            proteinSeq += map[codon2.replace('V', 'A')]
                        else:
                            proteinSeq += map[codon2.replace('V', 'C')]
                    elif 'S' in codon2:
                        if temp < 0.5:
                            proteinSeq += map[codon2.replace('S', 'G')]
                        else:
                            proteinSeq += map[codon2.replace('S', 'C')]
                    elif 'W' in codon2:
                        if temp < 0.5:
                            proteinSeq += map[codon2.replace('W', 'A')]
                        else:
                            proteinSeq += map[codon2.replace('W', 'T')]
                    elif 'D' in codon2:
                        if temp < 0.33:
                            proteinSeq += map[codon2.replace('D', 'G')]
                        elif temp > 0.67:
                            proteinSeq += map[codon2.replace('D', 'T')]
                        else:
                            proteinSeq += map[codon2.replace('D', 'A')]
                    elif 'Y' in codon2:
                        if temp < 0.5:
                            proteinSeq += map[codon2.replace('Y', 'T')]
                        else:
                            proteinSeq += map[codon2.replace('Y', 'C')]
                    elif 'R' in codon2:
                        if temp < 0.5:
                            proteinSeq += map[codon2.replace('R', 'C')]
                        else:
                            proteinSeq += map[codon2.replace('R', 'A')]
                    elif 'H' in codon2:
                        if temp < 0.33:
                            proteinSeq += map[codon2.replace('H', 'A')]
                        elif temp > 0.67:
                            proteinSeq += map[codon2.replace('H', 'T')]
                        else:
                            proteinSeq += map[codon2.replace('H', 'C')]

    ''' Print database query if applicable -- get result '''
    for item in proteins:
        if db.get(item):
            print db.get(item), '\n', '\t', item

db = collections.defaultdict(list)
def main():
    ''' Parse fasta file that is passed in -- parsing genome data '''
    f = fileinput.input()
    data = parseGenomeFile(f)

    ''' Setup database '''
    db = setupDB()

    ''' Translation'''
    for item in data.items():
        getProteinData(item[1])

if __name__ == "__main__":
    main()












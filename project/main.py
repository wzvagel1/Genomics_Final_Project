import fileinput
import collections
from collections import defaultdict
import os
import csv
import string

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
    for f in os.listdir('DB_Files/'):
        if f.endswith('.fa') or f.endswith('.fsa') or f.endswith('.fasta'):
            dbFiles.add(f)
    for f in os.listdir('.'):
        if f == 'db.csv':
            buildDB = False

    ''' Read in db.csv if it exists '''
    if not buildDB:
        first = True
        read = False
        for key, val in csv.reader(open("db.csv")):
            if first:
                if int(val) is len(dbFiles):
                    read = True
                first = False
            if read:
                db[key] = val
        else:
            if read:
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
        for key, val in db.items():
            w.writerow([key, val])
    return db

''' Taken from Dr. Langmead HW5 solutions '''
_revcomp_trans = string.maketrans("ACGTacgt", "TGCAtgca")
def getComplementStrand(s):
    return s[::-1].translate(_revcomp_trans)

def getAminoAcidSequence(strand):
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
                proteinSeq = []
                for i in xrange(x, len(seq), 3):
                    codon2 = seq[i:i+3]
                    if len(codon2) < 3: break
                    if codon2 in stop:
                        proteins.add(''.join(proteinSeq))
                        break
                    if codon2 in map.keys():
                        proteinSeq.append(map[codon2])

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
        proteins = getAminoAcidSequence(item[1])

if __name__ == "__main__":
    main()












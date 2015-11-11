import fileinput
import collections

map = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
    "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

def getProteins(strand, startPos):
    startCodon = "AUG"
    start = set()
    stop = set()
    temp = ""

    global proteins
    global invalidAlphabet
    for x in range(startPos, len(strand), 3):
        codon = strand[x:x+3]
        for c in codon:
            if c in invalidAlphabet:
                #codon contains invalid char
                print "CHECKING TEST"
                codon = ""
                break
        if len(codon) == 3:
            if map[codon] == "M":
                start.add(x)
            elif map[codon] == "STOP":
                stop.add(x)

    for item in sorted(start):
        toPrint = ""
        for x in range(item + 3, len(strand), 3):
            codon = strand[x-3: x]
            if map[codon] != "STOP":
                toPrint += map[codon]
            else:
                proteins.add(toPrint)
                break

proteins = set()
invalidAlphabet = {'B', 'D', 'E', 'F', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z'}
def main():
    f = fileinput.input()
    d = collections.defaultdict(list)
    temp = ""
    count = 0
    for i, line in enumerate(f):
        if '>' in line and i != 0:
            d[count] = temp
            temp = ""
            count += 1
        if '>' not in line: temp += line.rstrip()
    d[count] = temp
    ''' Modified so that multiple strands can be processed at once, just in case '''

    for s in d.items():
        strand = ""
        for char in s[1]:
            if char == "T":
                strand += "U"
            else: strand += char

        getProteins(strand, 0)
        getProteins(strand, 1)
        getProteins(strand, 2)

        ''' reading reverse complement '''
        complementStrand = ""
        for item in reversed(strand):
            if item == "C":
                complementStrand += 'G'
            if item == "G":
                complementStrand += 'C'
            if item == "A":
                complementStrand += 'U'
            if item == "U":
                complementStrand += 'A'

        getProteins(complementStrand, 0)
        getProteins(complementStrand, 1)
        getProteins(complementStrand, 2)

        print 'PROTEINS FOR ', s[0]
        for item in proteins:
            print '\t', item

if __name__ == "__main__":
    main()


















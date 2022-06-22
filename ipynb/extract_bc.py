import pyfastx
import sys

# python extract_bc.py fastq/LIB5456569_SAM24416644_S1_L001_R2_001.fastq.gz 25 80 wl/sense_wl_clipped.txt > output.txt


def read_wl(path, reads):
    unique = set()
    with open(path) as f:
        for bc in f:
            if bc.rstrip() in reads:
                unique.add(bc.rstrip())
    for i,bc in enumerate(unique):
        print(bc + ",BC" + str(i+1)) 


def parse_fastq(path, start, end):
    reads = set()
    fq = pyfastx.Fastx(path)
    for _,seq,__,___ in fq:
        reads.add(seq[start:end])
    return reads


if __name__ == "__main__":
    reads = parse_fastq(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
    read_wl(sys.argv[4], reads)


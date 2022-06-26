import pyfastx
import argparse


# python extract_bc.py --fastq fastq/LIB5456569_SAM24416644_S1_L001_R2_001.fastq.gz --start 25 --end 80 --whitelist wl/sense_wl_clipped.txt > intersected_bc.txt


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
    parser = argparse.ArgumentParser()
    parser.add_argument('--fastq', type=str, nargs=1, required=True,
                        help="Input read fastq with barcode segment")
    parser.add_argument('--start', type=int, nargs=1, required=True,
                        help="Start of the barcode segment, 0-indexed")
    parser.add_argument('--end', type=int, nargs=1,  required=True,
                        help="End of the barcode segment, 0-indexed")
    parser.add_argument('--whitelist', type=str, nargs=1, required=True,
                        help="Whitelist of all barcodes seperated by line")
    args = parser.parse_args()
    
    reads = parse_fastq(args.fastq, args.start, args.end)
    read_wl(args.whitelist, reads)

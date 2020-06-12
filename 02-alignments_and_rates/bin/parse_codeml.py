#!/usr/bin/python3
import os, sys, re

MODELS = ["null", "branch"]


def usage():
    print("usage: " + sys.argv[0] + " folder", file=sys.stderr)
    sys.exit(1)

def plausi():
    if len(sys.argv) != 2: usage()
    inFolder = sys.argv[1]
    return inFolder

def get_all_base_files(inFolder):
    fileHash = {}
    for file in os.listdir(inFolder):
        filename = os.path.split(file)[1]
        basename = filename.split("_")[0]
        fileHash[basename] = 1
    return fileHash.keys()

def parse_all_from_basefile(file, inFolder):
    modelHash = {}
    for m in MODELS:
        file_to_parse = inFolder + "/" + file + "_" + m + ".txt"

        if not os.path.exists(file_to_parse) or not os.path.isfile(file_to_parse):
            print("File not present: " + file_to_parse, file=sys.stderr)
            break

        # check if size of file is 0
        if os.stat(file_to_parse).st_size == 0:
            print('File is empty: ' + file_to_parse, file=sys.stderr)
            break

        fo = open(file_to_parse)
        for line in fo:
            if line.startswith("lnL("):
                np = re.match("lnL.*\s+np:\s*(\d+)", line ).group(1)
                lnL = re.match("lnL\(.*\):\s+([0-9.-]+)", line ).group(1)
                break
        try:
            lnL
        except NameError:
            print("File " + file_to_parse + " does not have the lnL field", file=sys.stderr)

        modelHash[ m ] = [lnL, np]
        fo.close()

    return modelHash

def main():
    inFolder  = plausi()
    basefiles = get_all_base_files(inFolder)
    all_results = {}
    for basefile in basefiles:
        all_results[ basefile ] = parse_all_from_basefile(basefile, inFolder)

    for gene in all_results:
        for model in all_results[ gene ]:
            sys.stdout.write(gene + "\t" + model + "\t" + all_results[ gene ][ model ][0] + "\t" + all_results[ gene ][ model ][1])
            sys.stdout.write("\n")

main()

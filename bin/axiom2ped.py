#!/usr/bin/env python3
import argparse
import csv
import collections
import logging
import os
import sys

from argparse import RawDescriptionHelpFormatter
from datetime import datetime
from collections import Counter
from shutil import copy2

################################# Argparse ###################################
parser = argparse.ArgumentParser(
    description="""Converts Axiom SNP Chip data and salmon individual information
        and outputs files in PED format.""",
    epilog =\
    f"Example of use: {sys.argv[0]}\
    -g SsaTrack_Island_1-8.calls_mod.txt\
    -a SsaTrack_Annotation_for_sharing.csv\
    -r SsaTrack_Island_riverInfo.csv\
    -i SsaTrack_Island_fishInfo.csv\
    -c Track_SNP_Pos.txt\
    -p Iceland\
    -o SsaTrack_Island_1-8", formatter_class=RawDescriptionHelpFormatter)

parser.add_argument("-g", "--genotypes", type=str, required=True,
                    help="Path to AXIOM file SsaTrack (e.g. SsaTrack_Island_1-8.calls_mod.txt")
parser.add_argument("-a", "--annotation", type=str, required=True,
                    help="Path to annotation file for SsaTrack (csv)")
parser.add_argument("-c", "--coordinates", type=str, required=True,
                    help="Probes coordinates (txt)")
parser.add_argument("-r", "--rivers", type=str, required=False,
                    help="Path to river information file (csv) (required for Icelandic Salmon)")                    
parser.add_argument("-i", "--individuals", type=str, required=False,
                    help="Path to individual information (csv) (required for Icelandic Salmon)")
parser.add_argument("-p", "--population", type=str, choices=['Iceland', 'Stofnfiskur'],
                    help = "Salmon population (Icelandic or Stofnfiskur", required = True)
parser.add_argument("-o", "--out_name", type=str, required=True,
                    help = "Name of output file")


if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()


if args.population == "Iceland":
    if args.individuals != None and args.rivers != None:
        pass
    else:
        sys.exit("ERROR: Since you are converting Icelandic samples, file containing information about rivers and individuals is required")


################################# Argparse ###################################

################################# Logging ####################################

logFileName = f"{parser.prog.strip('.py')}_{args.out_name}.log"
print(logFileName)
logfile_path = logFileName
# Remove existing log
if os.path.exists(logfile_path):
    os.remove(logfile_path)
else:
    pass

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

file_handler = logging.FileHandler(logfile_path)
file_handler.setFormatter(formatter)

logger.addHandler(file_handler)
#
logger.info(f"#{sys.argv[0]} was run with following argumetns")
logger.info(" ".join(sys.argv))
################################# Logging ####################################


def main():
    print("#INFO: Parsing genotype information")
    genotypes_dict = parse_ssatrack(args.genotypes)
    ssatrackAnnotation_dict = parse_ssatrackAnnotation(args.annotation)
    sstrackSNPcoordinates = parseCoordinates(args.coordinates)
    if args.population == "Iceland":
        print("#INFO: Population Iceland selected: Parsing individual information")
        individualsInfo_dict = parse_individualsInfo(args.individuals)
        riversInfo_dict = parse_riverInformation(args.rivers)
    elif args.population == "Stofnfiskur":
        print("#INFO: Population Stofnfiskur selected: No individual information presnent (expected)")
        individualsInfo_dict = parseStofnfiskurIndividuals(args.genotypes)
        riversInfo_dict = {}
    updateIdsFile = writePopFileForSNPdevelopment(individualsInfo_dict, args.out_name)
    write_plinkFiles(genotypes_dict, ssatrackAnnotation_dict, riversInfo_dict, individualsInfo_dict, args.out_name, args.population, sstrackSNPcoordinates)
    print(f"#INFO: To update population IDs to Ice, run plink with --update-ids {updateIdsFile}")
    print(f"Script exceuted successfully, check {logFileName}.log for more details")


def parse_ssatrack(genotypes):
    genotypes_dict = {}
    with open(genotypes, "r") as f_in:
        reader = csv.DictReader(filter(lambda row: row[0]!='#', f_in), delimiter='\t')
        for row in reader:
            probeSetId = row["probeset_id"]
            genotypes_dict[probeSetId] = row
    logger.info(f"Genotypic information of {len(genotypes_dict)} probes parsed from {genotypes}")
    return genotypes_dict


def parse_ssatrackAnnotation(annotation):
    ssatrackAnnotation_dict = {}
    with open(annotation, "r") as f_in:
        reader = csv.DictReader(filter(lambda row: row[0]!='#', f_in), delimiter=';')
        for row in reader:
            probeSetId = row['Probe_Set_ID']
            ssatrackAnnotation_dict[probeSetId] = row
    logger.info(f"{len(ssatrackAnnotation_dict)} annotations parsed from {annotation}")
    return(ssatrackAnnotation_dict)


def parseCoordinates(coordinates):
    sstrackSNPcoordinates  ={}
    no_positions = 0
    with open(coordinates, "r") as f_in:
        _ = f_in.readline()
        for line in f_in.readlines():
            sstrackSNPcoordinate = {}
            position = 0
            if "Mito" in line:
                probe, chromosome = line.strip().split()
                sstrackSNPcoordinate[probe] = {
                    "probe":probe,
                    "chromosome":"ssaMT",
                    "position":"0"
                }
                no_positions += 1
                sstrackSNPcoordinates[probe] = sstrackSNPcoordinate
            elif "unkn" in line:
                probe = line.strip().split()[0]
                sstrackSNPcoordinate[probe] = {
                    "probe":probe,
                    "chromosome":"0",
                    "position":"0"
                }
                no_positions += 1
                sstrackSNPcoordinates[probe] = sstrackSNPcoordinate
            else:
                probe, chromosome, position = line.strip().split()
                sstrackSNPcoordinate[probe] = {
                    "probe":probe,
                    "chromosome":chromosome,
                    "position":str(position)
                }
                no_positions += 1
                sstrackSNPcoordinates[probe] = sstrackSNPcoordinate
    logger.info(f"#Positions form {no_positions} parsed from {coordinates}")
    return(sstrackSNPcoordinates)
            


def parse_riverInformation(riversInfo):
    riversInfo_dict = {}
    with open(riversInfo, "r") as f_in:
        reader = csv.DictReader(f_in, delimiter=';')
        for row in reader:
            riverCode = row["river_code"]
            riversInfo_dict[riverCode] = row
    logger.info(f"Information about {len(riversInfo_dict)} rivers parsed from {riversInfo}")
    return(riversInfo_dict)


def parse_individualsInfo(individualsInfo):
    individualsInfo_dict = {}
    with open(individualsInfo, "r") as f_in:
        reader = csv.DictReader(f_in, delimiter=';')
        for row in reader:
            individualHash = row["Sample_hash"]
            individualsInfo_dict[individualHash] = row
    logger.info(f"Information about {len(individualsInfo_dict)} individuals parsed from {individualsInfo}")
    return(individualsInfo_dict)


def writePopFileForSNPdevelopment(individualsInfo_dict, out_name):
    # Function to write a txt file to change population ID of individuals (--update-ids)
    f_out_name = f"{out_name}_updateIds.txt"
    with open(f_out_name, "w") as f_out:
        for individual in individualsInfo_dict:
            # print(individual, individualsInfo_dict)
            old_pop = individualsInfo_dict[individual]["River_short"]
            individual_id = individual
            new_pop = "Ice"
            newPopInfoLine = f"{old_pop}\t{individual_id}\t{new_pop}\t{individual_id}\n"
            f_out.write(newPopInfoLine)
    return(f_out_name)


def parseStofnfiskurIndividuals(genotypes):
    individualsInfo_dict = {}
    individualsLine = ""
    with open(genotypes, "r") as f_in:
        for line in f_in.readlines():
            if line.split()[0] == "probeset_id":
                individualsLine = line
                break
        for individual in individualsLine.split()[1:]:
            individualsInfo_dict[individual.strip(".CEL")] = [0]
    logger.info(f"Information about {len(individualsInfo_dict)} individuals parsed from {genotypes}")
    return(individualsInfo_dict)


def write_plinkFiles(genotypes_dict, ssatrackAnnotation_dict, riversInfo_dict, individualsInfo_dict, out_name, population, sstrackSNPcoordinates):
    print("#INFO: Writing MAP File")
    probe_order = write_mapFile(ssatrackAnnotation_dict, genotypes_dict, out_name, sstrackSNPcoordinates)
    print("#INFO: Writing PED file")
    write_pedFile(genotypes_dict, ssatrackAnnotation_dict, riversInfo_dict, individualsInfo_dict, probe_order, out_name, population)
    return(probe_order)


def write_mapFile(ssatrackAnnotation_dict, genotypes_dict, out_name, sstrackSNPcoordinates):
    f_out_name = out_name + ".map"
    x = 0
    probesOrder = []
    with open(f_out_name, "w") as f_out:
        for probe in ssatrackAnnotation_dict:
            if probe in genotypes_dict.keys():
                chromosome = sstrackSNPcoordinates[probe][probe]["chromosome"]
                varID = probe
                pos_cm = "0"
                bp_coordinate = sstrackSNPcoordinates[probe][probe]["position"]
                f_out.write(f"{chromosome}\t{varID}\t{pos_cm}\t{bp_coordinate}\n")
                x += 1
                probesOrder.append(probe)
    logger.info(f"Map file written to {f_out_name}")
    return(probesOrder)


def write_pedFile(genotypes_dict, ssatrackAnnotation_dict, riversInfo_dict, individualsInfo_dict, probe_order, out_name, population):
    f_out_name = out_name + ".ped"
    with open(f_out_name, "w") as f_out:
        for individual in individualsInfo_dict:
            pedIndividInformation = generatepedIndividInformation(individual,riversInfo_dict,individualsInfo_dict, population) 
            pedGenotypeLine = []
            for probeSet in probe_order:
                axiom_genotype = genotypes_dict[probeSet][individual + ".CEL"]
                ped_genotype = axiom2ped_genotype(axiom_genotype, genotypes_dict, probeSet, ssatrackAnnotation_dict, individual)
                pedGenotypeLine.append(ped_genotype)
            pedGenotypeString = '\t'.join(pedGenotypeLine)
            ped_line = f"{pedIndividInformation}\t{pedGenotypeString}\n"
            f_out.write(ped_line)
    logger.info(f"Ped file written to {f_out_name}")


def generatepedIndividInformation(individual, riversInfo_dict, individualsInfo_dict, population):
    pedIndividInformation = ""
    if population == "Iceland":
        river_code = individualsInfo_dict[individual]["River_short"]
        population = river_code
        pedSex = parseSex(individualsInfo_dict[individual]["Sex"])
    elif population == "Stofnfiskur":
        population = "Stofn"
        pedSex = "0"
    paternalID, maternalID, sex, phenotype = "0", "0", pedSex, "-9"
    pedIndividInformation = "\t".join([population, individual, paternalID, maternalID, sex, phenotype])
    return(pedIndividInformation)


def parseSex(sex_string):
    ped_sex_dict = {"M":"1", "F":"2", "0":"0"}
    return(ped_sex_dict[sex_string])


def axiom2ped_genotype(axiom_genotype, genotypes_dict, probeSet, ssatrackAnnotation_dict, individual):
    ped_genotype = ""
    alleleA, alleleB = ssatrackAnnotation_dict[probeSet]["Allele_A"], ssatrackAnnotation_dict[probeSet]["Allele_B"]
    axiomGenotype = genotypes_dict[probeSet][individual + ".CEL"]
    axiom2ped_dict = {
        "-1":"0 0",
        "0":f"{alleleA} {alleleA}",
        "1":f"{alleleA} {alleleB}",
        "2":f"{alleleB} {alleleB}"
    }
    ped_genotype = axiom2ped_dict[axiomGenotype]
    return(ped_genotype)


### NEED ALSO TO PARSE STOFNFISKUR, SHOULD IT BE IN THE SAME OR SEPERATE SCRIPT.
### PROBABLY BEST TO HAVE A FLAG IN ARGUMENTS, WITH EITHER ICELAND OR BOTH

main()

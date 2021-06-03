#!/usr/bin/env python3
import argparse
import csv
import logging
import os
import sys

from argparse import RawDescriptionHelpFormatter
from shutil import copy2

################################# Argparse ###################################
parser = argparse.ArgumentParser(
    description="""Converts Axiom SNP Chip data and salmon individual information
        and outputs files in PED format.""",
    epilog =\
    f"Example of use: {sys.argv[0]}\
    -s SsaTrack_Island_1-8.calls_mod.txt\
    -f SsaTrack_Sten_Island.calls_mod.txt\
    -a SsaTrack_Annotation_for_sharing.csv\
    -r SsaTrack_Island_riverInfo.csv\
    -i SsaTrack_Island_fishInfo.csv\
    -o SsaTrack_Island_Stofnfiskur_1-8", formatter_class=RawDescriptionHelpFormatter)

parser.add_argument("-s", "--ssatrack", type=str, required=True,
                    help="Path to Icelandic AXIOM file SsaTrack")
parser.add_argument("-f", "--stofnfiskur", type=str, required=False,
                    help="Path to STOFNFISKUR AXIOM file SsaTrack")
parser.add_argument("-a", "--annotation", type=str, required=True,
                    help="Path to annotation file for SsaTrack (csv)")
parser.add_argument("-r", "--rivers", type=str, required=True,
                    help="Path to river information file (csv)")                    
parser.add_argument("-i", "--individuals", type=str, required=True,
                    help="Path to individual information (csv)")
parser.add_argument("-o", "--out_name", type=str, required=True,
                    help = "Name of output file")

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()


def check_if_file_exists():
    for inputFile in [args.ssatrack, args.stofnfiskur, args.annotation, args.rivers, args.individuals]:
        if not os.path.exists(inputFile):
            parser.error(f"The file {inputFile} does not exist!")
        

check_if_file_exists()

################################# Argparse ###################################

################################# Logging ####################################

logfile_path = f"axiom2ped.log"
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
    if args.stofnfiskur is None:
        print("#INFO: Parsing only Icelandic samples")
    ssatrack_dict = parse_ssatrack(args.ssatrack)
    ssatrackAnnotation_dict = parse_ssatrackAnnotation(args.annotation)
    riversInfo_dict = parse_riverInformation(args.rivers)
    individualsInfo_dict = parse_individualsInfo(args.individuals)
    write_plinkFiles(ssatrack_dict, ssatrackAnnotation_dict, riversInfo_dict, individualsInfo_dict, args.out_name)
    if args.stofnfiskur is not None:
        print("#INFO Parsing both Icelandic and Stofnfiskur samples")
        ssatrack_stofnfiskur_dict = parse_stofnfiskur_ssatrack(args.stofnfiskur)
        stofniskur_individuals = parse_stofnfiskur_inividuals(args.stofnfiskur)
        addStofnfiskurPlinkFiles(ssatrack_stofnfiskur_dict, stofniskur_individuals, args.out_name, ssatrackAnnotation_dict)
    print("Script exceuted successfully, check axiom2ped.log for more details")


def parse_ssatrack(ssatrack):
    ssatrack_dict = {}
    with open(ssatrack, "r") as f_in:
        reader = csv.DictReader(filter(lambda row: row[0]!='#', f_in), delimiter='\t')
        for row in reader:
            probeSetId = row["probeset_id"]
            ssatrack_dict[probeSetId] = row
    logger.info(f"Genotypic information of {len(ssatrack_dict)} probes parsed from {ssatrack}")
    return ssatrack_dict


def parse_stofnfiskur_ssatrack(stofnfiskur):
    stofnfiskur_dict = {}
    with open(stofnfiskur, "r") as f_in:
        reader = csv.DictReader(filter(lambda row: row[0]!='#', f_in), delimiter='\t')
        for row in reader:
            probeSetId = row["probeset_id"]
            stofnfiskur_dict[probeSetId] = row
    logger.info(f"Genotypic information of {len(stofnfiskur_dict)} probes parsed from {stofnfiskur}")
    return stofnfiskur_dict


def parse_stofnfiskur_inividuals(stofnfiskur):
    stofniskur_individuals = []
    with open(stofnfiskur, "r") as f_in:
        for line in f_in.readlines():
            if line.split("\t")[0] == "probeset_id":
                for individ in line.strip().split("\t"):
                    if individ != "probeset_id":
                        stofniskur_individuals.append(individ.replace(".CEL",""))
    logger.info(f"{len(stofniskur_individuals)} parsed from {stofnfiskur}")
    return stofniskur_individuals


def parse_ssatrackAnnotation(annotation):
    ssatrackAnnotation_dict = {}
    with open(annotation, "r") as f_in:
        reader = csv.DictReader(filter(lambda row: row[0]!='#', f_in), delimiter=';')
        for row in reader:
            probeSetId = row['Probe_Set_ID']
            ssatrackAnnotation_dict[probeSetId] = row
    logger.info(f"{len(ssatrackAnnotation_dict)} annotations parsed from {annotation}")
    return(ssatrackAnnotation_dict)


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


def write_plinkFiles(ssatrack_dict, ssatrackAnnotation_dict, riversInfo_dict, individualsInfo_dict, out_name):
    probe_order = write_mapFile(ssatrackAnnotation_dict, ssatrack_dict, out_name)
    write_pedFile(ssatrack_dict, ssatrackAnnotation_dict, riversInfo_dict, individualsInfo_dict, probe_order, out_name)
    return(probe_order)


def write_mapFile(ssatrackAnnotation_dict, ssatrack_dict, out_name):
    f_out_name = out_name + ".map"
    x = 0
    probesOrder = []
    with open(f_out_name, "w") as f_out:
        for probe in ssatrackAnnotation_dict:
            if probe in ssatrack_dict.keys():
                probe_info = ssatrackAnnotation_dict[probe]
                chromosome = ""
                if ssatrackAnnotation_dict[probe]["Chromosome"] == "MT":
                    chromosome = "MT"
                else:
                    chromosome = "0"
                varID = probe
                pos_cm, bp_coordinate = "0", "0"
                f_out.write(f"{chromosome}\t{varID}\t{pos_cm}\t{bp_coordinate}\n")
                x += 1
                probesOrder.append(probe)
    logger.info(f"Map file written to {f_out_name}")
    return(probesOrder)


def write_pedFile(ssatrack_dict, ssatrackAnnotation_dict, riversInfo_dict, individualsInfo_dict, probe_order, out_name):
    f_out_name = out_name + ".ped"
    with open(f_out_name, "w") as f_out:
        for individual in individualsInfo_dict:
            pedIndividInformation = generatepedIndividInformation(individual,riversInfo_dict,individualsInfo_dict) 
            pedGenotypeLine = []
            for probeSet in probe_order:
                axiom_genotype = ssatrack_dict[probeSet][individual + ".CEL"]
                ped_genotype = axiom2ped_genotype(axiom_genotype, ssatrack_dict, probeSet, ssatrackAnnotation_dict, individual)
                pedGenotypeLine.append(ped_genotype)
            pedGenotypeString = '\t'.join(pedGenotypeLine)
            ped_line = f"{pedIndividInformation}\t{pedGenotypeString}\n"
            f_out.write(ped_line)
    logger.info(f"Ped file written to {f_out_name}")


def generatepedIndividInformation(individual, riversInfo_dict, individualsInfo_dict):
    pedIndividInformation = ""
    river_code = individualsInfo_dict[individual]["River_short"]
    population = river_code
    pedSex = parseSex(individualsInfo_dict[individual]["Sex"])
    paternalID, maternalID, sex, phenotype = "0", "0", pedSex, "-9"
    pedIndividInformation = "\t".join([population, individual, paternalID, maternalID, sex, phenotype])
    return(pedIndividInformation)

def parseSex(sex_string):
    ped_sex_dict = {"M":"1", "F":"2", "0":"0"}
    return(ped_sex_dict[sex_string])


def axiom2ped_genotype(axiom_genotype, ssatrack_dict, probeSet, ssatrackAnnotation_dict, individual):
    ped_genotype = ""
    alleleA, alleleB = ssatrackAnnotation_dict[probeSet]["Allele_A"], ssatrackAnnotation_dict[probeSet]["Allele_B"]
    axiomGenotype = ssatrack_dict[probeSet][individual + ".CEL"]
    axiom2ped_dict = {
        "-1":"0 0",
        "0":f"{alleleA} {alleleA}",
        "1":f"{alleleA} {alleleB}",
        "2":f"{alleleB} {alleleB}"
    }
    ped_genotype = axiom2ped_dict[axiomGenotype]
    return(ped_genotype)


def addStofnfiskurPlinkFiles(ssatrack_stofnfiskur_dict, stofniskur_individuals, out_name, ssatrackAnnotation_dict):
    probeOrderMapFile = getProbeOrderMapFile(out_name)
    addToPedFile(ssatrack_stofnfiskur_dict, stofniskur_individuals, out_name, probeOrderMapFile, ssatrackAnnotation_dict)


def getProbeOrderMapFile(out_name):
    probeOrder = []
    mapFileIn = out_name + ".map"
    with open(mapFileIn, "r") as f_in:
        for line in f_in.readlines():
            _, probe, _, _ = line.strip().split("\t")
            probeOrder.append(probe)
    return(probeOrder)
            

def addToPedFile(ssatrack_stofnfiskur_dict, stofniskur_individuals, out_name, probeOrderMapFile, ssatrackAnnotation_dict):
    f_ped_in_out = out_name + ".ped"
    with open(f_ped_in_out, "a") as f_out:
        for stofnfisk in stofniskur_individuals:
            pedIndividInformation = f"Stofn\t{stofnfisk}\t0\t0\t0\t-9"
            pedGenotypeLine = []
            for probeSet in probeOrderMapFile:
                axiom_genotype = ssatrack_stofnfiskur_dict[probeSet][stofnfisk + ".CEL"]
                ped_genotype = axiom2ped_genotype(axiom_genotype, ssatrack_stofnfiskur_dict, probeSet, ssatrackAnnotation_dict, stofnfisk)
                pedGenotypeLine.append(ped_genotype)
            pedGenotypeString = '\t'.join(pedGenotypeLine)
            ped_line = f"{pedIndividInformation}\t{pedGenotypeString}\n"
            f_out.write(ped_line)
    logger.info(f"Stofnfiskur samples added to existing PLINK file")
    print("#Script finished without errors")




main()

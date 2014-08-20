from __future__ import division
import sys
import os
from collections import defaultdict
import gzip

###################################################################################################

def zopen(fileName, mode):
    if fileName.endswith(".gz"):
        return gzip.open(fileName, mode)
    else:
        return open(fileName, mode)

###################################################################################################

def checkVCF(vcfName):

    os.system("vcf-validator %s" %(vcfName))
    os.system("python scripts/vcfChecks.py %s" %(vcfName))

    if vcfName.endswith("gz"):
        os.system("zcat %s | grep -v REFCALL | python scripts/computeTsTv.py" %(vcfName))
        os.system("zcat %s | grep -v REFCALL | python scripts/binIndelRatioByHP.py 3" %(vcfName))
    else:
        os.system("cat %s | grep -v REFCALL | python scripts/computeTsTv.py" %(vcfName))
        os.system("cat %s | grep -v REFCALL | python scripts/binIndelRatioByHP.py 3" %(vcfName))

    vcfFile = zopen(vcfName, 'r')
    filters = defaultdict(int)
    varTypes = defaultdict(int)
    passVarTypes = defaultdict(int)
    genotypes = defaultdict(int)
    inconsistentGenotypes = defaultdict(int)

    for line in vcfFile:
        if line.startswith("#"):
            continue

        cols = line.strip().split("\t")
        ref = cols[3]
        alts = cols[4]
        theFilters = cols[6]
        sampleInfo = cols[9]

        for theFilter in theFilters.split(";"):
            filters[theFilter] += 1

        for alt in alts.split(","):

            if alt == ".":
                continue

            if len(ref) == len(alt):

                if len(ref) == 1:
                    varTypes["SNP"] += 1
                else:
                    varTypes["MNP"] += 1

            elif len(ref) > len(alt):
                varTypes['Deletion'] += 1
                varTypes['Indel'] += 1
            else:
                varTypes['Insertion'] += 1
                varTypes['Indel'] += 1

            if theFilters == "PASS":
                if len(ref) == len(alt):

                    if len(ref) == 1:
                        passVarTypes["SNP"] += 1
                    else:
                        passVarTypes["MNP"] += 1

                elif len(ref) > len(alt):
                    passVarTypes['Deletion'] += 1
                    passVarTypes['Indel'] += 1
                else:
                    passVarTypes['Insertion'] += 1
                    passVarTypes['Indel'] += 1

        sampleInfo = cols[9].split(":")
        genotype = sampleInfo[0]
        genotypes[genotype] += 1

        # Only look at likelihoods etc for bi-allelic sites
        if not "," in alts and alts != ".":
            nv = int(sampleInfo[-1])
            nr = int(sampleInfo[-2])
            likelihoods = sampleInfo[1].split(",")
            maxLikeIndex = likelihoods.index(max(likelihoods))

            if genotype == "0/0":
                if maxLikeIndex != 0:
                    inconsistentGenotypes["Hom Ref LL/GT Mismatch"] += 1
                if nv > 0:
                    inconsistentGenotypes["Hom Ref With >0 Var Reads"] += 1

                if maxLikeIndex == 0 and nv == 0:
                    inconsistentGenotypes["Hom Ref Ok"] += 1

            elif genotype == "0/1" or genotype == "1/0":

                if maxLikeIndex != 1:
                    inconsistentGenotypes["Het LL/GT Mismatch"] += 1

                if nv/nr < 0.15:
                    inconsistentGenotypes["Het With < 15% VarReads"] += 1

                if nv/nr > 0.85:
                    inconsistentGenotypes["Het With > 85% VarReads"] += 1

                if maxLikeIndex == 1 and 0.15 <= nv/nr <= 0.95:
                    inconsistentGenotypes["Het Ok"] += 1

            elif genotype == "1/1":

                if maxLikeIndex != 2:
                    inconsistentGenotypes["Hom Var LL/GT Mismatch"] += 1

                if nv/nr < 0.9:
                    inconsistentGenotypes["Hom Var With < 90% VarReads"] += 1

                if maxLikeIndex == 2 and nv/nr > 0.9:
                    inconsistentGenotypes["Hom Var Ok"] += 1


    vcfFile.close()
    totalGenotypes = sum(genotypes.values())
    nInconsistentGenotypes = sum(inconsistentGenotypes.values())

    print ""
    print "Filters:"

    for theFilter,count in filters.iteritems():
        print "%s : %s" %(theFilter, count)

    print ""
    print "varTypes:"

    for theType,count in varTypes.iteritems():
        print "%s : %s" %(theType, count)

    print ""
    print "passVarTypes:"

    for theType,count in passVarTypes.iteritems():
        print "%s : %s" %(theType, count)

    print ""
    print "Genotype Summary:"
    for key,val in genotypes.iteritems():
        print "%s = %s (%1.2f %%)" %(key,val,100*val/totalGenotypes)

    print ""
    print "Genotype Consistency Summary:"
    for key,val in inconsistentGenotypes.iteritems():
        print "%s = %s (%1.2f %%)" %(key,val,100*val/nInconsistentGenotypes)

    #print ""
    #print "1Kg Membership:"

    #if vcfName.endswith("gz"):
    #    os.system("zcat %s | python scripts/computePhaseOneMembership.py Phase1SNPs_Chr20.gz" %(vcfName))
    #else:
    #    os.system("cat %s | python scripts/computePhaseOneMembership.py Phase1SNPs_Chr20.gz" %(vcfName))

###################################################################################################

if len(sys.argv) < 5:
    print ""
    print "Invalid usage."
    print ""
    print "Correct usage as follows:"
    print "python sanityChecks.py PATH/Platypus.py data.bam ref.fa output.vcf OTHER_ARGS"
    print ""
    print ""
    sys.exit(1)

platypus = sys.argv[1]
bams = sys.argv[2]
ref = sys.argv[3]
vcfName = sys.argv[4]
rest = " ".join(sys.argv[5:])

command = "time python %s callVariants --bamFiles=%s --refFile=%s --output=%s --regions=20 %s" %(platypus, bams,ref,vcfName,rest)
os.system(command)
checkVCF(vcfName)

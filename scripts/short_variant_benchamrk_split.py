import sys, os, subprocess
# from dataclasses import dataclass
from collections import defaultdict
# import pandas as pd
# import numpy as np


dir_name = sys.argv[1]
threshold = sys.argv[2]

# @dataclass
# class Gene:
#     tp: int = 0
#     fp: int = 0
#     fn: int = 0
amp = 1
max_quality = "QUAL>"+threshold
min_quality = "QUAL<"+threshold
first = 1
# df = pd.DataFrame(columns = ["Gene", "TP", "STP", "ITP", "ITPi", "ITPd",
#                                      "EFN", "SEFN", "IEFN", "IEFNi", "IEFNd",
#                                      "FP", "SFP", "IFP", "IFPi", "IFPd",
#                                      "FN", "SFN", "IFN", "IFNi", "IFNd",])
header = ["Gene", "TP", "STP", "ITP", "ITPi", "ITPd",
                                     "EFN", "SEFN", "IEFN", "IEFNi", "IEFNd",
                                     "FP", "SFP", "IFP", "IFPi", "IFPd",
                                     "FN", "SFN", "IFN", "IFNi", "IFNd",
                                     "APR","ARE","AF1","SPR","SRE","SF1","IPR","IRE","IF1","IPRI","IREI","IF1I","IPRD","IRED","IF1D"]
print("\t".join(header))

# Note: for indels one indel could be fiting into insertion and deletions e.g., chr19	51361697	.	GTT	G,GTTT
for subdir, dirs, files in os.walk(dir_name):
    if first:
        first = 0
        continue
    gene_name = (os.path.basename(os.path.normpath(subdir))).split('__')[0]
    genes_benchmark = defaultdict(list)
    for file in files:
        if file in ["tp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz"]:
            if file == "tp.vcf.gz":
                tp = int(subprocess.run(["bcftools view -H -i '" + max_quality + "' " + os.path.join(subdir, file) + " | wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n")) # We have '' around the filter input otherwise it will not work.
                stp = int(subprocess.run(["bcftools view -H -v snps -i '" + max_quality + "' " + os.path.join(subdir, file) + " | wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n"))
                itp = int(subprocess.run(["bcftools view -v indels -i '" + max_quality + "' " + os.path.join(subdir, file) + " | bcftools norm -m - | bcftools view -H | wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n"))
                itpi = int(subprocess.run(["bcftools view -v indels -i '" + max_quality + "' "+ os.path.join(subdir, file) + "| bcftools norm -m - |  bcftools filter --include 'strlen(REF)<strlen(ALT)' |   bcftools view -H | wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n"))
                itpd = int(subprocess.run(["bcftools view  -v indels -i '" + max_quality + "' " + os.path.join(subdir, file) + " | bcftools norm -m - |   bcftools filter --include 'strlen(REF)>strlen(ALT)' | bcftools view -H| wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n"))
                genes_benchmark[gene_name].append(tp)
                genes_benchmark[gene_name].append(stp)
                genes_benchmark[gene_name].append(itp)
                genes_benchmark[gene_name].append(itpi)
                genes_benchmark[gene_name].append(itpd)
                efn = int(subprocess.run(["bcftools view -H -i '" + min_quality + "' " + os.path.join(subdir, file) + " | wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n")) # We have '' around the filter input otherwise it will not work.
                sefn = int(subprocess.run(["bcftools view -H -v snps -i '" + min_quality + "' " + os.path.join(subdir, file) + " | wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n"))
                iefn = int(subprocess.run(["bcftools view -v indels -i '" + min_quality + "' " + os.path.join(subdir, file) + " | bcftools norm -m - | bcftools view -H | wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n"))
                iefni = int(subprocess.run(["bcftools view -v indels -i '" + min_quality + "' "+ os.path.join(subdir, file) + "| bcftools norm -m - |  bcftools filter --include 'strlen(REF)<strlen(ALT)' |   bcftools view -H | wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n"))
                iefnd = int(subprocess.run(["bcftools view  -v indels -i '" + min_quality + "' " + os.path.join(subdir, file) + " | bcftools norm -m - |   bcftools filter --include 'strlen(REF)>strlen(ALT)' | bcftools view -H| wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n"))
                genes_benchmark[gene_name].append(efn)
                genes_benchmark[gene_name].append(sefn)
                genes_benchmark[gene_name].append(iefn)
                genes_benchmark[gene_name].append(iefni)
                genes_benchmark[gene_name].append(iefnd)
            elif file == "fp.vcf.gz":
                fp = int(subprocess.run(["bcftools view -H -i '" + max_quality + "' " + os.path.join(subdir, file) + " | wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n"))
                sfp = int(subprocess.run(["bcftools view -H -v snps -i '" + max_quality + "' " + os.path.join(subdir, file) + " | wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n"))
                ifp = int(subprocess.run(["bcftools view -v indels -i '" + max_quality + "' " + os.path.join(subdir, file) + " | bcftools view -H | wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n"))
                ifpi = int(subprocess.run(["bcftools view -v indels -i '" + max_quality + "' "+ os.path.join(subdir, file) + " |  bcftools filter --include 'strlen(REF)<strlen(ALT)' | bcftools view -H | wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n"))
                ifpd = int(subprocess.run(["bcftools view  -v indels -i '" + max_quality + "' " + os.path.join(subdir, file) + " | bcftools filter --include 'strlen(REF)>strlen(ALT)' | bcftools view -H | wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n"))
                genes_benchmark[gene_name].append(fp)
                genes_benchmark[gene_name].append(sfp)
                genes_benchmark[gene_name].append(ifp)
                genes_benchmark[gene_name].append(ifpi)
                genes_benchmark[gene_name].append(ifpd)
            elif file == "fn.vcf.gz":
                fn = int(subprocess.run(["bcftools view -H " + os.path.join(subdir, file) + " | wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n"))
                sfn = int(subprocess.run(["bcftools view -H -v snps -i '" + max_quality + "' " + os.path.join(subdir, file) + " | wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n"))
                ifn = int(subprocess.run(["bcftools view -v indels -i '" + max_quality + "' " + os.path.join(subdir, file) + " | bcftools view -H | wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n"))
                ifni = int(subprocess.run(["bcftools view -v indels -i '" + max_quality + "' "+ os.path.join(subdir, file) + " |  bcftools filter --include 'strlen(REF)<strlen(ALT)' | bcftools view -H | wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n"))
                ifnd = int(subprocess.run(["bcftools view  -v indels -i '" + max_quality + "' " + os.path.join(subdir, file) + " | bcftools filter --include 'strlen(REF)>strlen(ALT)' | bcftools view -H | wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n"))
                genes_benchmark[gene_name].append(fn)
                genes_benchmark[gene_name].append(sfn)
                genes_benchmark[gene_name].append(ifn)
                genes_benchmark[gene_name].append(ifni)
                genes_benchmark[gene_name].append(ifnd)
            if len(genes_benchmark[gene_name]) == 20:#len(header) - 1:
                # adding false negative resuls based on theshold we used
                fn = fn + efn
                sfn = sfn + sefn
                ifn = ifn + iefn
                ifni = ifni + iefni
                ifnd = ifnd + iefnd

                # calculating precison, recall and F1-score:
                # apr=tp/(tp+fp)*amp  # Precision

                apr=(tp+fp)*amp and tp/(tp+fp)*amp or -1 # Precision
                # are=tp/(tp+fn)*amp # Recall # Note: it eother the recall is really zero or there is no variants
                are=(tp+fn)*amp and tp/(tp+fn)*amp or -1  # a / b # Recall # Note: it eother the recall is really zero or there is no variants
                # af1=2*apr*are/(apr+are) # F1-score
                af1=(apr+are) and 2*apr*are/(apr+are) or -1  # a / b # F1-score
                # spr=stp/(stp+sfp)*amp # Substitution Precision
                spr=(stp+sfp)*amp and stp/(stp+sfp)*amp or -1 # Substitution Precision
                # sre=stp/(stp+sfn)*amp # Substitution Recall
                sre= (stp+sfn)*amp and stp/(stp+sfn)*amp or -1  # a / b# Substitution Recall
                # sf1=2*spr*sre/(spr+sre) # Substitution F1-score
                sf1= (spr+sre) and 2*spr*sre/(spr+sre) or -1 # Substitution F1-score
                ipr=ire=if1=0

                # Initiate variable in case if no indel true positive itp.
                ipr=-1
                ire=-1
                if1=-1
                ipri=-1
                irei=-1
                if1i=-1
                iprd=-1
                ired=-1
                if1d=-1
                
                if itp > 0:
                  # ipr=itp/(itp+ifp)*amp
                  # ire=itp/(itp+ifn)*amp
                  # if1=2*ipr*ire/(ipr+ire)
                  # ipri=itpi/(itpi+ifpi)*amp
                  # irei=itpi/(itpi+ifni)*amp
                  ipr=(itp+ifp)*amp and itp/(itp+ifp)*amp or -1
                  ire=(itp+ifn)*amp and itp/(itp+ifn)*amp or -1
                  if1=(ipr+ire) and 2*ipr*ire/(ipr+ire) or -1
                  ipri=(itpi+ifpi)*amp and itpi/(itpi+ifpi)*amp or -1
                  irei=(itpi+ifni)*amp and itpi/(itpi+ifni)*amp or -1
                  if1i=(ipri+irei) and 2*ipri*irei/(ipri+irei) or -1  # a / b
                  # if1i=2*ipri*irei/(ipri+irei)  (ipri+irei) and 2*ipri*irei/(ipri+irei) or -1  # a / b
                  # iprd=itpd/(itpd+ifpd)*amp
                  iprd= (itpd+ifpd)*amp and itpd/(itpd+ifpd)*amp or -1
                  # ired=itpd/(itpd+ifnd)*amp
                  ired=(itpd+ifnd)*amp and itpd/(itpd+ifnd)*amp or -1
                  # if1d=2*iprd*ired/(iprd+ired)
                  if1d=(iprd+ired) and 2*iprd*ired/(iprd+ired) or -1


                # [apr, are, af1, spr, sre, sf1, ipr, ire, if1, ipri, irei, if1i, iprd, ired, if1d]
                s = "\t".join([str(i) for i in genes_benchmark[gene_name]])
                print(f"{gene_name}\t{s}\t{apr}\t{are}\t{af1}\t{spr}\t{sre}\t{sf1}\t{ipr}\t{ire}\t{if1}\t{ipri}\t{irei}\t{if1i}\t{iprd}\t{ired}\t{if1d}")
            # if len(genes_benchmark[gene_name]) == len(df.columns) -1:
            #     df.loc[len(df.index)] = [gene_name, *genes_benchmark[gene_name]]

# print(df)

# We can use normalization like:
# ifp = subprocess.run(["bcftools view -v indels -i '" + max_quality + "' " + os.path.join(subdir, file) + " | bcftools norm -m - | bcftools view -H | wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n")
# But we will not to avoid complexity
#---------------------------------------

# Devide by zero withput error
# result = b and a / b or 0  # a / b
# When b != 0 we have True and a / b or 0. True and a / b is equal to a / b. a / b or 0 is equal to a / b.
# When b == 0 we have False and a / b or 0. False and a / b is equal to False. False or 0 is equal to 0.

# for k, v in genes_benchmark.items():
#     print(f"{k} {v}")


#
# print("-------------------")
# for filename in glob.iglob(dir_name+"/**/[t|f][p|n].vcf.gz", recursive=True):
#     if os.path.isfile(filename): # filter dirs
#         print(filename)

# if file in ["tp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz"]:
#     if gene_name in genes_benchmark:
#         pass
#     else:
#         genes_benchmark[gene_name] = Gene(0, 0, 0)
#
#     if file == "tp.vcf.gz":
#         tp = subprocess.run(["bcftools view -H -i '" + max_quality + "' " + os.path.join(subdir, file) + " | wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n") # We have '' around the filter input otherwise it will not work.
#         genes_benchmark[gene_name].tp = tp
#     elif file == "fp.vcf.gz":
#         fp = subprocess.run(["bcftools view -H -i '" + max_quality + "' " + os.path.join(subdir, file) + " | wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n")
#         genes_benchmark[gene_name].fp = fp
#     elif file == "fn.vcf.gz":
#         fn = subprocess.run(["bcftools view -H " + os.path.join(subdir, file) + " | wc -l" ], capture_output=True, encoding='UTF-8', shell=True, text=True).stdout.strip("\n")
#         genes_benchmark[gene_name].fn = fn

#  1	Gene Gene name
#  2	TP True positive
#  3	STP Substitution True poistive
#  4	ITP Indels True positive`
#  5	ITPi Indels true positive insertions
#  6	ITPd Indels true positive deletions
#  7	EFN  filter value False negative
#  8	SEFN filter value Substitution False negative
#  9	IEFN
# 10	IEFNi
# 11	IEFNd
# 12	FP
# 13	SFP
# 14	IFP
# 15	IFPi
# 16	IFPd
# 17	FN
# 18	SFN
# 19	IFN
# 20	IFNi
# 21	IFNd
# 22	APR    # all Precision
# 23	ARE    # all Recall
# 24	AF1    # all F1-score
# 25	SPR    # Substitution
# 26	SRE
# 27	SF1
# 28	IPR    # Indels
# 29	IRE
# 30	IF1
# 31	IPRI
# 32	IREI
# 33	IF1I
# 34	IPRD
# 35	IRED
# 36	IF1D

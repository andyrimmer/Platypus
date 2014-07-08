#
# add allele count field
#

import vcf
import filez
import pysam
import sys
import getopt
import math
import bisect


min_matching_stretch = 999    # require a dimer at minimum

tandem_thresholds =      [0,6,9, 11,13,14,16,18,18,18]                             # about 1% of loci polymorphic
tandem_thresholds_high = [0,11,15,22,27,30,31,31,31,31,31,31,31,31,31,31,31,31]    # about 30% of loci polymorphic


# set to True if collecting background palindromic counts; otherwise collects background indel rates
do_palindromes = False

# set to True if creating a mask for highly mutable regions; otherwise prints records
# set to True if getting G/C region lengths
do_mask = True
do_get_gcregion = False

# masked regions will be merged if they are closer than mask_merge nucleotides
mask_merge = 10

# coefficients for indel rate model from YRI calls
_mu =  [0.40, 0.08, 0.019,0.016,0.01, 0.011,0.01, 0.014,0.011,0.015,0.015,0.012,0.01,0.013,0.008,0.01, 0.006,0.006,0.004,0.006,0.003,0.006,0.005,0.004,0.03]
_r =  [0.13,0.33, 3.4, 6.1,9.4, 37, 134,98, 172,328,470,286,297,179,104,479,276,5,80,80,80,80,80,80,80,80,80,80]

# indel error model; from AW_SC_4654.bam 
indel_error_model={'AAAAAAGAC': 'VVVVVVVVVVVVVCCC', 1: 'VVVMJGB<83/,*)))))))))))))))))))))))))))', 2: 'VVVNNMLIGB>:8653311/.-,++++**)))))))))))(((((((((((((((((', 3: 'VVVVVMLKJHFCA?><<;9976654433330000000/.................((((((((((((', 4: 'VVVVVVVKKJHEECA@??<;::9876664443322100//////////..------+++*(((((((', 5: 'VVVVVVVVVKJIGEBAA?><;;;;9888885533221000----------------,,,,,', 6: 'VVVVVVVVVVJJHHFDBAA?>>=<<:::8876632221111111-', 7: 'VVVVVVVVVVVIIHGECBA@?>=<;;;;;76', 8: 'VVVVVVVVVVVVIEEEDCBA?><<<<<::::77777777777777777777777777777777777****))', 9: 'VVVVVVVVVVVVVIEEDDCAA?==<<88888774', 10: 'VVVVVVVVVVVVVVHFECAA@@><;;997777773', 11: 'VVVVVVVVVVVVVVVHGEDCBA?>>;;:::7', 'AGC': 'VVVVVNNNMLJDDC?<<;98', 'AAACTT': 'VVVVVVVVVVMMH', 'ACCCC': 'VVVVVVVVVGEEE', 'AAATT': 'VVVVVVVVVLKJI', 'AAACTG': 'VVVVVVVVVVKK', 'ACATGT': 'VVVVVVVVVVH', 'AAAATTT': 'VVVVVVVVVVVIGGGGC', 'AAAAAAGAG': 'VVVVVVVVVVVVVCCCB', 'ACACATACAT': 'VVVVVVVVVVVVVVVV>', 'AAACCTGAC': 'VVVVVVVVVVVVVVVH', 'AAAAATAAAT': 'VVVVVVVVVVVVVVHEBBBA=', 'AATCAG': 'VVVVVVVVVVMMI', 'AAAAATAT': 'VVVVVVVVVVVVDDA', 'AAAACAT': 'VVVVVVVVVVVJJI', 'AAAAAAATAC': 'VVVVVVVVVVVVVV>====', 'AAAAAAAGAAG': 'VVVVVVVVVVVVVVV@', 'AACAGG': 'VVVVVVVVVVK', 'AAAAAGAAAG': 'VVVVVVVVVVVVVVVB===::', 'AAAAATAC': 'VVVVVVVVVVVVDAAA?', 'AAAACAG': 'VVVVVVVVVVVJ', 'AAAAATAG': 'VVVVVVVVVVVVIE', 'AAATATTT': 'VVVVVVVVVVVVH', 'AAAAAAATAT': 'VVVVVVVVVVVVVVA>', 'ATCCCCC': 'VVVVVVVVVVVF', 'CCCCG': 'VVVVVVVVVIGA', 'AGCGGGG': 'VVVVVVVVVVVVB', 'AAATCC': 'VVVVVVVVVVI', 'AAATCT': 'VVVVVVVVVVI', 'AAAGTC': 'VVVVVVVVVVL', 'AGCAT': 'VVVVVVVVVK', 'ACACCAGC': 'VVVVVVVVVVVVJI', 'ACTGCGCC': 'VVVVVVVVVVVVD', 'AAAAAAAAT': 'VVVVVVVVVVVVV=====<<<<<8', 'AAAATAAT': 'VVVVVVVVVVVVHFFFB?', 'ACAGAG': 'VVVVVVVVVVIICCCCC>', 'AACAG': 'VVVVVVVVVKKI', 'AATAGCT': 'VVVVVVVVVVVE', 'AACAC': 'VVVVVVVVVKKK', 'AACTAG': 'VVVVVVVVVVJJ', 'AAAACT': 'VVVVVVVVVVHHD', 'AAAATAAG': 'VVVVVVVVVVVVHHH', 'ACAGAT': 'VVVVVVVVVVJI', 'AAAGACC': 'VVVVVVVVVVVJ', 'AAAACC': 'VVVVVVVVVVGG', 'AACAT': 'VVVVVVVVVJJJ', 'AACCTGGCC': 'VVVVVVVVVVVVVI', 'AAACCC': 'VVVVVVVVVVJ', 'AAGGAGGGAG': 'VVVVVVVVVVVVVVV@', 'AGATCC': 'VVVVVVVVVVJJ', 'ACATACATAT': 'VVVVVVVVVVVVVVD', 'AGATCG': 'VVVVVVVVVVVF', 'AAAAACAAAAC': 'VVVVVVVVVVVVVVVBBBBBAAAA:', 'AAACCT': 'VVVVVVVVVVG', 'AGAGGGG': 'VVVVVVVVVVVGG', 'AATCTCTG': 'VVVVVVVVVVVVLL', 'AAAAACC': 'VVVVVVVVVVVGGF', 'AGCATCC': 'VVVVVVVVVVVG', 'ACCCATC': 'VVVVVVVVVVVG', 'AAAAATAAC': 'VVVVVVVVVVVVVVB', 'AGGG': 'VVVVVVVIHFFDC>>>;;99985553', 'AAAAATAAG': 'VVVVVVVVVVVVVG', 'AGGC': 'VVVVVVVKKJJJH', 'AAAAAAATTTT': 'VVVVVVVVVVVVVVVB?', 'AAAATAATAAT': 'VVVVVVVVVVVVVVV@', 'AAAAAAGAAG': 'VVVVVVVVVVVVVVCCB', 'AGCCTCCTC': 'VVVVVVVVVVVVVH', 'AAAATAAAGT': 'VVVVVVVVVVVVVVVVVV?', 'ACCTGG': 'VVVVVVVVVVI', 'AAAATGTAT': 'VVVVVVVVVVVVVG', 'ACCTGC': 'VVVVVVVVVVH', 'AAAAATAAT': 'VVVVVVVVVVVVVGFFE', 'AAATGAT': 'VVVVVVVVVVVL', 'AAGGGAGG': 'VVVVVVVVVVVVEDC', 'AAAGGGGG': 'VVVVVVVVVVVVD', 'ACATATGTAT': 'VVVVVVVVVVVVVVJ', 'AATGT': 'VVVVVVVVVJHHH', 'ACACCC': 'VVVVVVVVVVG', 'ACACACACC': 'VVVVVVVVVVVVVVVV>', 'AAATAGT': 'VVVVVVVVVVVH', 'AATGG': 'VVVVVVVVVIGE@=;;99777755411111000+))))))$$$$$$$$$$$$$$$$$$$', 'AATGC': 'VVVVVVVVVP', 'AAAGCAG': 'VVVVVVVVVVVI', 'AAAGAC': 'VVVVVVVVVVNMI', 'AAAGAG': 'VVVVVVVVVVLKKJJJA?', 'AAAAAAAAAG': 'VVVVVVVVVVVVVV766666666664', 'AAAAATTAC': 'VVVVVVVVVVVVVVBBB', 'AGAGAGG': 'VVVVVVVVVVVH', 'AAAGAT': 'VVVVVVVVVVKG', 'AACAGT': 'VVVVVVVVVVK', 'AAATTGG': 'VVVVVVVVVVVL', 'AAAGAGAGAG': 'VVVVVVVVVVVVVVDA@@@@<', 'AAAAAAAGAG': 'VVVVVVVVVVVVVVB>', 'AATGTGCAC': 'VVVVVVVVVVVVVVVM', 'AAGGAAGGAGG': 'VVVVVVVVVVVVVVVVV>', 'ACAGAGTG': 'VVVVVVVVVVVVVH', 'ACATCC': 'VVVVVVVVVVI', 'ACCATC': 'VVVVVVVVVVHFFB>:8887', 'AAAAAGCT': 'VVVVVVVVVVVVFF', 'AAAAAGAAC': 'VVVVVVVVVVVVVG', 'AAAAGTC': 'VVVVVVVVVVVI', 'AAAAGTG': 'VVVVVVVVVVVI', 'AAGGAGGG': 'VVVVVVVVVVVVEE?===', 'AAAAGTT': 'VVVVVVVVVVVGF', 'AGCTCC': 'VVVVVVVVVVL', 'AAAAAGAAAAG': 'VVVVVVVVVVVVVVV@?>=<;;;;7', 'AAGAGG': 'VVVVVVVVVVKKKH', 'ACATC': 'VVVVVVVVVL', 'AAGAGC': 'VVVVVVVVVVI', 'AAAGAAAG': 'VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV,++++))', 'AAATGC': 'VVVVVVVVVVM', 'AAACAAT': 'VVVVVVVVVVVJJH', 'AAATGG': 'VVVVVVVVVVHHHF', 'ACCCT': 'VVVVVVVVVJI', 'AAAAAACAAAC': 'VVVVVVVVVVSVVVVVBBAAAA????7', 'AAACAAG': 'VVVVVVVVVVVJ', 'AAAGAAAT': 'VVVVVVVVVVVVIIH', 'AAATGT': 'VVVVVVVVVVIIIG', 'AATCAC': 'VVVVVVVVVVL', 'AAAAAAATTT': 'VVVVVVVVVVVVVVBBBBA', 'AAATGAC': 'VVVVVVVVVVVJ', 'AAAGTG': 'VVVVVVVVVVLLJ', 'AAAAAATC': 'VVVVVVVVVVVVF', 'AAAAAATG': 'VVVVVVVVVVVVEEEC', 'ACCATGCGCC': 'VVVVVVVVVVVVVVE', 'ACACCTC': 'VVVVVVVVVVVF', 'AAAGTT': 'VVVVVVVVVVJ', 'AAAAAATT': 'VVVVVVVVVVVVFFEDA', 'AAAAGGT': 'VVVVVVVVVVVF', 'AAGTGG': 'VVVVVVVVVVK', 'AAGTGT': 'VVVVVVVVVVI', 'AAACAGC': 'VVVVVVVVVVVG', 'AGATGG': 'VVVVVVVVVVJ', 'AAAATGGGAG': 'VVVVVVVVVVVVVVI', 'AAAGAATG': 'VVVVVVVVVVVVE', 'AGCC': 'VVVVVVVLLKJ', 'AAAAAAGAAAG': 'VVVVVVVVVVVVVVVCA??<<::9', 'AGCG': 'VVVVVVVL', 'AACTCC': 'VVVVVVVVVVG', 'ACCACCATGC': 'VVVVVVVVVVVVVVHG', 'ACCTCC': 'VVVVVVVVVVJJH', 'AACT': 'VVVVVVVJJIFE', 'AGCT': 'VVVVVVVOOKF', 'A': 'VVVMJGC<83/,**))))))))))))))))))))))))))', 'AATATT': 'VVVVVVVVVVI', 'AAGTC': 'VVVVVVVVVL', 'ACGAGATC': 'VVVVVVVVVVVVB', 'AAATACC': 'VVVVVVVVVVVF', 'AATCT': 'VVVVVVVVVJJ', 'AATATG': 'VVVVVVVVVVJ', 'AATCC': 'VVVVVVVVVLK', 'AACATGGC': 'VVVVVVVVVVVVJI', 'AAATACT': 'VVVVVVVVVVVG', 'AAGGAC': 'VVVVVVVVVVJ', 'AAGGAGGGAGG': 'VVVVVVVVVVVVVVVVV?', 'AAGGAG': 'VVVVVVVVVVKH', 'AAAAAACC': 'VVVVVVVVVVVVEEB?', 'ACACTGGG': 'VVVVVVVVVVVVA', 'AAATGAAT': 'VVVVVVVVVVVVJ', 'AGCCTGGCC': 'VVVVVVVVVVVVVH', 'AAAAAAACC': 'VVVVVVVVVVVVV@@@', 'AAGGAT': 'VVVVVVVVVVI', 'AAATGAAC': 'VVVVVVVVVVVVE', 'AAAAAACT': 'VVVVVVVVVVVVF', 'ACAG': 'VVVVVVVMKKKJECCA?', 'AAAAGAT': 'VVVVVVVVVVVHGG', 'AAAT': 'VVVVVVVKIIHFFDBB??>=;;:9766554333210//////////////////', 'AAAATTC': 'VVVVVVVVVVVJ', 'AAGAGAGG': 'VVVVVVVVVVVVFFF', 'AAAACAAAC': 'VVVVVVVVVVVVVDDCA', 'AAAC': 'VVVVVVVKKKFDC@@@@?<;:998666653333222222', 'AAGCCTC': 'VVVVVVVVVVVF', 'ACCTCT': 'VVVVVVVVVVI', 'AAAG': "VVVVVVVKKIGAA@?<<;:777755330//....****))))(((((('''''''''''''''''''", 'ACCCTC': 'VVVVVVVVVVH', 'AGATC': 'VVVVVVVVVL', 'AAAGGAT': 'VVVVVVVVVVVGG', 'AAAAATG': 'VVVVVVVVVVVJH', 'AAAGGAAG': 'VVVVVVVVVVVVHHD', 'ACCCAGCC': 'VVVVVVVVVVVVC', 'AAAAATC': 'VVVVVVVVVVVJ', 'AAAGGAAC': 'VVVVVVVVVVVVF', 'AATATGT': 'VVVVVVVVVVVVE', 'AATGGGG': 'VVVVVVVVVVVH', 'AAAGGAG': 'VVVVVVVVVVVMKI', 'AAAAATT': 'VVVVVVVVVVVJFFDC', 'ACCCAGCT': 'VVVVVVVVVVVVF', 'ACATAGAT': 'VVVVVVVVVVVVVV>', 'AAAAAATTTTT': 'VVVVVVVVVVVVVVVCCCB', 'AACAAG': 'VVVVVVVVVVJ', 'ACACAG': 'VVVVVVVVVVL', 'AAGAG': 'VVVVVVVVVMKKKG', 'AAAATACAC': 'VVVVVVVVVVVVVC', 'AAAAGAGG': 'VVVVVVVVVVVVH', 'ACACAT': 'VVVVVVVVVVIHDDDDB', 'AACAAT': 'VVVVVVVVVVI', 'AATGTG': 'VVVVVVVVVVL', 'AAGAT': 'VVVVVVVVVKJJ', 'ACACACACGC': 'VVVVVVVVVVVVVVVVVV;', 'AAAATACAG': 'VVVVVVVVVVVVVD', 'AAAATCT': 'VVVVVVVVVVVIF', 'AAAACTAC': 'VVVVVVVVVVVVJJ', 'AAAAAACAC': 'VVVVVVVVVVVVVAAAA', 'AAGGTGG': 'VVVVVVVVVVVH', 'AAAGGTGGGG': 'VVVVVVVVVVVVVVVVA', 'AAAACTAT': 'VVVVVVVVVVVVGGGG', 'AAAATCC': 'VVVVVVVVVVVHHH', 'ACTC': 'VVVVVVVOLJJEEE<', 'AAAATTG': 'VVVVVVVVVVVI', 'ACTG': 'VVVVVVVNNJJJE', 'AAGTCT': 'VVVVVVVVVVJ', 'AAAAAAAGG': 'VVVVVVVVVVVVVV@?', 'ACAGAGCGAG': 'VVVVVVVVVVVVVVF', 'AGAGGGAGGG': 'VVVVVVVVVVVVVV@@@', 'ACACACATGC': 'VVVVVVVVVVVVVVBBBA>', 'ACTCAGCT': 'VVVVVVVVVVVVM', 'AACTATCTAC': 'VVVVVVVVVVVVVVK', 'ACACCCC': 'VVVVVVVVVVVE', 'AAAAAAATGG': 'VVVVVVVVVVVVVV>>', 'AACATAT': 'VVVVVVVVVVVG', 'AAAAAATAC': 'VVVVVVVVVVVVVB>>>>>', 'AATGAT': 'VVVVVVVVVVKE', 'AAAAAATAG': 'VVVVVVVVVVVVVEAAA', 'AAAGGAAGG': 'VVVVVVVVVVVVVEEEEE:', 'AAGTG': 'VVVVVVVVVLKJ', 'AATGAC': 'VVVVVVVVVVK', 'AAAAAATAT': 'VVVVVVVVVVVVVC@', 'AATGAG': 'VVVVVVVVVVM', 'AACTGGCT': 'VVVVVVVVVVVVVD', 'ACCTCCC': 'VVVVVVVVVVVJIF', 'AAAAAAACAAC': 'VVVVVVVVVVVVVVV?', 'AAACAGAC': 'VVVVVVVVVVVVGC', 'AAAATAAC': 'VVVVVVVVVVVVH', 'CCGG': 'VVVVVVVG', 'AAGGTG': 'VVVVVVVVVVLK', 'AAAAAGGG': 'VVVVVVVVVVVVE', 'AAGGGAGGG': 'VVVVVVVVVVVVVC', 'ACATATATAT': 'VVVVVVVVVVVVVVC??<<;;;;;;;;66', 'AGCCCC': 'VVVVVVVVVVHH', 'AAAAAAGG': 'VVVVVVVVVVVVEB', 'AATGGATG': 'VVVVVVVVVVVVVVED', 'ACTGGAG': 'VVVVVVVVVVVG', 'AAAACAACAAC': 'VVVVVVVVVVVVVVVC', 'AAAAAAGT': 'VVVVVVVVVVVVFFDB', 'ATCCC': 'VVVVVVVVVJII', 'ACTCTCT': 'VVVVVVVVVVVG', 'AAAAAATTTT': 'VVVVVVVVVVVVVVGGGD?', 'ACCTCCCTC': 'VVVVVVVVVVVVV?', 'ACCACT': 'VVVVVVVVVVI', 'AAAAAAAGTT': 'VVVVVVVVVVVVVVV>', 'ATC': 'VVVVVONMLJJEB???<<:9887754444431111111-', 'ACATGCC': 'VVVVVVVVVVVF', 'AAACTGG': 'VVVVVVVVVVVG', 'AAAATTAT': 'VVVVVVVVVVVVH', 'AATAATG': 'VVVVVVVVVVVJ', 'AAATAAGT': 'VVVVVVVVVVVVVCCC', 'AATAATC': 'VVVVVVVVVVVK', 'AAAAAATGT': 'VVVVVVVVVVVVVD', 'ACTCCCC': 'VVVVVVVVVVVHE', 'AATGGTGC': 'VVVVVVVVVVVVVVH', 'AATAATT': 'VVVVVVVVVVVG', 'AAAGGT': 'VVVVVVVVVVJ', 'AAATCAG': 'VVVVVVVVVVVI', 'AAAAGTAG': 'VVVVVVVVVVVVF', 'AACTGC': 'VVVVVVVVVVK', 'AAAAACAT': 'VVVVVVVVVVVVEC', 'AAAACAAAAT': 'VVVVVVVVVVVVVVEEEEE', 'AAAGGG': 'VVVVVVVVVVJH', 'AAAGGC': 'VVVVVVVVVVI', 'AAAAACT': 'VVVVVVVVVVVJF', 'AAAAACAC': 'VVVVVVVVVVVVCCBB', 'AAAAATAATTT': 'VVVVVVVVVVVVVVVD', 'AACCCCC': 'VVVVVVVVVVVE', 'AAAAATTTTT': 'VVVVVVVVVVVVVVCCB', 'AAAAT': 'VVVVVVVVVHHHFFCAAA@@==<<<977775552.........---------', 'AAAATGT': 'VVVVVVVVVVVJH', 'AATC': 'VVVVVVVMLLLKKKKK@>>>>>>8884', 'AATCTAT': 'VVVVVVVVVVVG', 'AATG': 'VVVVVVVNNNKHHEDCBAAA;:::::7766655551', 'AAATTTC': 'VVVVVVVVVVVK', 'AAAATGC': 'VVVVVVVVVVVHF', 'AAAAC': 'VVVVVVVVVIIHFFCCBAA@====::888877443332110', 'AATCAGGC': 'VVVVVVVVVVVVII', 'AAATGAATG': 'VVVVVVVVVVVVVVVVVVI', 'AAAATGG': 'VVVVVVVVVVVJH', 'AAAAG': 'VVVVVVVVVIHHE@===<;77777777555555555,', 'AAG': 'VVVVVNMLKJGA?>>;:88865444443333333333333333333333333333((((((((((((', 'AAC': 'VVVVVMMLIHFBB@==<:::9666654433000+++++****&', 'AAAAGAAAAG': 'VVVVVVVVVVVVVVVVV7', 'ACCG': 'VVVVVVVG', 'AAATAATAAT': 'VVVVVVVVVVVVVVVB', 'AGGGGGC': 'VVVVVVVVVVVHFD', 'AAAGAAAGG': 'VVVVVVVVVVVVVBB?;;;', 'AAGACTG': 'VVVVVVVVVVVG', 'AACTGG': 'VVVVVVVVVVJ', 'AGGGCC': 'VVVVVVVVVVII', 'AAAGCAGG': 'VVVVVVVVVVVVII', 'ACCT': 'VVVVVVVLJJ', 'AAAGCAGC': 'VVVVVVVVVVVVIH', 'AACAGGGC': 'VVVVVVVVVVVVVC', 'AAAATATTT': 'VVVVVVVVVVVVVIIGD', 'ACTCCT': 'VVVVVVVVVVI', 'AAGGAGAG': 'VVVVVVVVVVVVD', 'AACAGAGC': 'VVVVVVVVVVVVVE', 'AAAGTAG': 'VVVVVVVVVVVJI', 'ACTCCG': 'VVVVVVVVVV@', 'ACTCCC': 'VVVVVVVVVVI', 'AAGCAG': 'VVVVVVVVVVO', 'AAGCAC': 'VVVVVVVVVVI', 'AAAAACAAC': 'VVVVVVVVVVVVVFDB@', 'AAAGATTT': 'VVVVVVVVVVVVHG', 'AAATGGGCT': 'VVVVVVVVVVVVVVH', 'AAGCAT': 'VVVVVVVVVVJ', 'CCCG': 'VVVVVVVHH', 'AAAAATATT': 'VVVVVVVVVVVVVHHG', 'AAATAAATT': 'VVVVVVVVVVVVVFFE', 'ACAGGAG': 'VVVVVVVVVVVH', 'AACCC': 'VVVVVVVVVK', 'ACAGGC': 'VVVVVVVVVVK', 'AAAAAT': 'VVVVVVVVVVHHEECCCCCA?????>', 'AAAGAAC': 'VVVVVVVVVVVI', 'AAAGAAG': 'VVVVVVVVVVVKGD', 'ACTCTGC': 'VVVVVVVVVVVE', 'ACCTCCCC': 'VVVVVVVVVVVVVV@', 'AAAAAC': 'VVVVVVVVVVHHGFDDDBA@>>======999995', 'AAAATACC': 'VVVVVVVVVVVVE', 'AAAAAAACAC': 'VVVVVVVVVVVVVV@>', 'AACCT': 'VVVVVVVVVLF', 'AAAAAG': 'VVVVVVVVVVHFFDDA?<<<;;;7', 'AAAATACG': 'VVVVVVVVVVVVB', 'AGGAGGCCG': 'VVVVVVVVVVVVVGG', 'AAAGAAT': 'VVVVVVVVVVVK', 'ACATCCTCTCC': 'VVVVVVVVVVVVVVVVH', 'AAACAAAG': 'VVVVVVVVVVVVHHG', 'AAAAAAAATAC': 'VVVVVVVVVVVVVVV77777', 'AACAATAAT': 'VVVVVVVVVVVVVC', 'ACCAGT': 'VVVVVVVVVVLI', 'AAAATGTTT': 'VVVVVVVVVVVVVF', 'AACAGGT': 'VVVVVVVVVVVD', 'AAACAAAT': 'VVVVVVVVVVVVHHD', 'AGATAT': 'VVVVVVVVVVVVVVVV???????????.......+', 'ACCTGAGATC': 'VVVVVVVVVVVVVVBB', 'AGGGGGGC': 'VVVVVVVVVVVVA@', 'AATTAT': 'VVVVVVVVVVIG', 'AAAACAAAAG': 'VVVVVVVVVVVVVVD', 'AAACCAG': 'VVVVVVVVVVVG', 'AAAAGAAT': 'VVVVVVVVVVVVJ', 'AATTAC': 'VVVVVVVVVVKF', 'AATTAG': 'VVVVVVVVVVLLJ', 'AAAAGAG': 'VVVVVVVVVVVH', 'ACTCTG': 'VVVVVVVVVVH', 'AAAAGAAG': 'VVVVVVVVVVVVFD', 'AATGCCTGT': 'VVVVVVVVVVVVVF', 'AAATTC': 'VVVVVVVVVVOK', 'AAATTG': 'VVVVVVVVVVL', 'AAATACAT': 'VVVVVVVVVVVVG', 'ACACATGC': 'VVVVVVVVVVVVVKID', 'AATAGAATGG': 'VVVVVVVVVVVVVVVVV<:', 'ACCTC': 'VVVVVVVVVJ', 'ACCTG': 'VVVVVVVVVIII', 'AAATTT': 'VVVVVVVVVVJI', 'ACCACCATC': 'VVVVVVVVVVVVVBBB999888222/+', 'AGAGGAGG': 'VVVVVVVVVVVVE', 'ACCTCCCCC': 'VVVVVVVVVVVVV>', 'AACATG': 'VVVVVVVVVVL', 'AATGGCGC': 'VVVVVVVVVVVVI', 'AAAAGGAAAG': 'VVVVVVVVVVVVVVDAA?<:', 'AAATG': 'VVVVVVVVVNLHHHDA', 'AAAATT': 'VVVVVVVVVVJIIIIECB', 'AAAATTTT': 'VVVVVVVVVVVVIG', 'AAATGAGG': 'VVVVVVVVVVVVF', 'ACTGG': 'VVVVVVVVVML', 'ACCTCCTG': 'VVVVVVVVVVVVVVC', 'AAAATACAT': 'VVVVVVVVVVVVVF', 'ACTGC': 'VVVVVVVVVKKJJ', 'ACCCCCC': 'VVVVVVVVVVVECC', 'AAAATG': 'VVVVVVVVVVIIG', 'AAAATTTC': 'VVVVVVVVVVVVI', 'AAAATC': 'VVVVVVVVVVKIIH', 'AAAATTTG': 'VVVVVVVVVVVVI', 'AAGG': 'VVVVVVVLKJHFDB<<<<<99999755533332211+***)))))))))))))))))))((', 'AAGC': 'VVVVVVVLJJ', 'ACTCCACTGC': 'VVVVVVVVVVVVVVE', 'AAATGCAG': 'VVVVVVVVVVVVM', 'AAGAGGG': 'VVVVVVVVVVVI', 'AAAGAAAGAG': 'VVVVVVVVVVVVVV???==::::', 'AGGGGCC': 'VVVVVVVVVVVH', 'ACGC': 'VVVVVVVIIG', 'AAGT': 'VVVVVVVJII', 'AAATTAT': 'VVVVVVVVVVVJH', 'ATCCCC': 'VVVVVVVVVVHHH', 'ACGG': 'VVVVVVVIB', 'AGGGGC': 'VVVVVVVVVVK', 'AAAGGAGG': 'VVVVVVVVVVVVF', 'AAATCAAT': 'VVVVVVVVVVVVHGF', 'AGGGGG': 'VVVVVVVVVVHFBBA', 'AAATAAATAAT': 'VVVVVVVVVVVVVVVVVVVVVVV;', 'AAAAATCT': 'VVVVVVVVVVVVD', 'AGATG': 'VVVVVVVVVOONNN?????5', 'AAAATGGC': 'VVVVVVVVVVVVFF', 'ATGCCCC': 'VVVVVVVVVVVE', 'AACATT': 'VVVVVVVVVVK', 'AAAACCT': 'VVVVVVVVVVVF', 'AAAAAGAT': 'VVVVVVVVVVVVH', 'AAAAATCC': 'VVVVVVVVVVVVEEC', 'AAGAAGG': 'VVVVVVVVVVVJ', 'ACCAT': 'VVVVVVVVVMKG', 'AAAAAGAC': 'VVVVVVVVVVVVH', 'ACACTC': 'VVVVVVVVVVVVVVVV;', 'ACCCTG': 'VVVVVVVVVVI', 'AAAAAGAG': 'VVVVVVVVVVVVIIH', 'AAAACCC': 'VVVVVVVVVVVVC', 'AAAACCCC': 'VVVVVVVVVVVVA', 'C': 'VVVLIE?84-,,,,,,+', 'AAATAC': 'VVVVVVVVVVIG', 'AAATAG': 'VVVVVVVVVVNK', 'ATGCC': 'VVVVVVVVVJ', 'AGCCC': 'VVVVVVVVVKKJH', 'AAAGAAGG': 'VVVVVVVVVVVVFCAAAAAAAA6', 'AAATAT': 'VVVVVVVVVVKHF', 'AATT': 'VVVVVVVLLLIIIGGGGGGG==<', 'AGGCG': 'VVVVVVVVVV?', 'ACCAG': 'VVVVVVVVVLL', 'AGCCT': 'VVVVVVVVVI', 'ACACACACAT': 'VVVVVVVVVVVVVVAAAA@@?', 'AAATGGG': 'VVVVVVVVVVVI', 'AGGCGG': 'VVVVVVVVVVJHH', 'ACAGCC': 'VVVVVVVVVVM', 'AAAAAAATT': 'VVVVVVVVVVVVV????', 'ACCCCT': 'VVVVVVVVVVH', 'AAACAG': 'VVVVVVVVVVNNIH', 'AAACAC': 'VVVVVVVVVVJ', 'ACCCCC': 'VVVVVVVVVVFF', 'AAAACAAT': 'VVVVVVVVVVVVGGE', 'AAACAT': 'VVVVVVVVVVIIIG', 'AAAAAAG': 'VVVVVVVVVVVEEDDB?=======5', 'AAAAGTGT': 'VVVVVVVVVVVVHH', 'AAAGAGGTAC': 'VVVVVVVVVVVVVVI', 'AAAAAAC': 'VVVVVVVVVVVFFFC@@@@>>>>>>::', 'AAAACTT': 'VVVVVVVVVVVG', 'AAAAATACT': 'VVVVVVVVVVVVVDDA?=', 'AGCTC': 'VVVVVVVVVO', 'AGGGGGG': 'VVVVVVVVVVVCCCC@>>', 'AAAAAATAAAT': 'VVVVVVVVVVVVVVVEC?', 'AAAGAAGTAG': 'VVVVVVVVVVVVVVC', 'AAAGGGAG': 'VVVVVVVVVVVVG', 'AAAAATACC': 'VVVVVVVVVVVVVVA@', 'AAAAATACG': 'VVVVVVVVVVVVVV?', 'AATACT': 'VVVVVVVVVVH', 'ACACC': 'VVVVVVVVVHGGGEE5', 'AAAAGAAAT': 'VVVVVVVVVVVVVIIG', 'CG': 'VVVNLKHFC?;9954', 'ACCCACTGC': 'VVVVVVVVVVVVVC', 'ATATC': 'VVVVVVVVVN', 'ACACT': 'VVVVVVVVVL', 'AAAAAATGAT': 'VVVVVVVVVVVVVVC', 'ACCC': 'VVVVVVVJIGFDAA<', 'AAAAGAAAG': 'VVVVVVVVVVVVVF@>>==<:', 'AAATAATT': 'VVVVVVVVVVVVIIF', 'AAGGGG': 'VVVVVVVVVVJH', 'AAAACAAG': 'VVVVVVVVVVVVH', 'ACTAT': 'VVVVVVVVVII', 'AAT': 'VVVVVLKIIHFD@@><<;997644332211000000////..........-', 'AAGGAGGCC': 'VVVVVVVVVVVVVGF', 'AAGAAGGAAGG': 'VVVVVVVVVVVVVVV?', 'ACTCT': 'VVVVVVVVVK', 'AGAGGC': 'VVVVVVVVVVLI', 'ACCAGCCTG': 'VVVVVVVVVVVVVII', 'AGAGGG': 'VVVVVVVVVVJGGGG@99999999985', 'AAAAAAAAG': 'VVVVVVVVVVVVV<<;;;;;;:', 'AAAACAAC': 'VVVVVVVVVVVVVC', 'ACTCC': 'VVVVVVVVVHHHHCCCCC2--', 'AACTAAT': 'VVVVVVVVVVVG', 'AAAAACAAAC': 'VVVVVVVVVVVVVVBBB??', 'ACTGCT': 'VVVVVVVVVVH', 'AACC': 'VVVVVVVLII', 'ACCCCTG': 'VVVVVVVVVVVEEEE??', 'AAATTAG': 'VVVVVVVVVVVG', 'ACCCCTC': 'VVVVVVVVVVVG', 'ACTGCC': 'VVVVVVVVVVKJ', 'ACACACCACC': 'VVVVVVVVVVVVVVD', 'ACACACAT': 'VVVVVVVVVVVVGEEED', 'AACCACTTG': 'VVVVVVVVVVVVVE', 'AATATAT': 'VVVVVVVVVVVHC', 'AAAAAATTC': 'VVVVVVVVVVVVVVVBB', 'AAAGAGAG': 'VVVVVVVVVVVVGDAAA>=', 'AAAGAAAGAAG': 'VVVVVVVVVVVVVVVBB??::::8', 'AAAGAGAC': 'VVVVVVVVVVVVF', 'AAAAAATTT': 'VVVVVVVVVVVVVDD', 'AACTT': 'VVVVVVVVVK', 'ATGC': 'VVVVVVVNN', 'AATACAT': 'VVVVVVVVVVVH', 'ACACACACAG': 'VVVVVVVVVVVVVVAAAA', 'AAGCT': 'VVVVVVVVVK', 'AACTG': 'VVVVVVVVVP', 'AACTC': 'VVVVVVVVVNM', 'ATATATATATC': 'VVVVVVVVVVVVVVVVV;;;;;5', 'AAGGAGG': 'VVVVVVVVVVVJJI', 'AATAATAT': 'VVVVVVVVVVVVH', 'AGGCCC': 'VVVVVVVVVVJ', 'AAAATAG': 'VVVVVVVVVVVH', 'AGGCCG': 'VVVVVVVVVVJIG', 'AAAATAC': 'VVVVVVVVVVVAA', 'AAGGGGG': 'VVVVVVVVVVVF', 'AAAAATTTT': 'VVVVVVVVVVVVVGGEEB', 'ACCATGCC': 'VVVVVVVVVVVVFFF', 'AGG': 'VVVVVLLJJIFE?>><;::86663111/.', 'AACGTGCAC': 'VVVVVVVVVVVVVL', 'AAAATAT': 'VVVVVVVVVVVII', 'AAGAC': 'VVVVVVVVVL', 'AAAAAAGTC': 'VVVVVVVVVVVVVG', 'AAATGCT': 'VVVVVVVVVVVH', 'AGAGAGGG': 'VVVVVVVVVVVVHDDD', 'ACACATAT': 'VVVVVVVVVVVVEB@@', 'AAGGCAG': 'VVVVVVVVVVVJ', 'ACACACACACC': 'VVVVVVVVVVVVVVVVVVV8', 'AAAAAGAAG': 'VVVVVVVVVVVVVHD', 'ACAGCCC': 'VVVVVVVVVVVG', 'AACTTG': 'VVVVVVVVVVKK', 'CCG': 'VVVVVIIIIBBBB======:::1', 'AAATTTCC': 'VVVVVVVVVVVVG', 'AAGTAT': 'VVVVVVVVVVI', 'AGGCC': 'VVVVVVVVVIG', 'AAAAATAAAAT': 'VVVVVVVVVVVVVVVFDAAA=', 'AAAGTGG': 'VVVVVVVVVVVGG', 'ACATAT': 'VVVVVVVVVVIGDDDCCCCC9', 'AAAATGTT': 'VVVVVVVVVVVVH', 'AAAAAATCC': 'VVVVVVVVVVVVVV?', 'ACAGAGAG': 'VVVVVVVVVVVVH', 'AAAAATATTTT': 'VVVVVVVVVVVVVVVD', 'AAAATGTG': 'VVVVVVVVVVVVH', 'AAAAAATACC': 'VVVVVVVVVVVVVVV@', 'AGCCTC': 'VVVVVVVVVVKII', 'CCCCCG': 'VVVVVVVVVVDDD', 'ACCTCAGGTG': 'VVVVVVVVVVVVVVE', 'AAGTGATCTC': 'VVVVVVVVVVVVVVC', 'AAAAATGC': 'VVVVVVVVVVVVFB', 'AAAATCTAG': 'VVVVVVVVVVVVVVH', 'AAAAATGG': 'VVVVVVVVVVVVH', 'AAATAAT': 'VVVVVVVVVVVIIE', 'ACAGC': 'VVVVVVVVVPO', 'ACAGG': 'VVVVVVVVVL', 'AAGTGATTC': 'VVVVVVVVVVVVVKH', 'ACCACTGC': 'VVVVVVVVVVVVJ', 'AAATAAG': 'VVVVVVVVVVVIIG', 'AAAAATGT': 'VVVVVVVVVVVVHCC', 'ACAGT': 'VVVVVVVVVL', 'AATCTG': 'VVVVVVVVVVM', 'AAAAAAAC': 'VVVVVVVVVVVV@@@@@?=<<;;;:::::::+', 'AATAC': 'VVVVVVVVVL', 'AATCTC': 'VVVVVVVVVVJ', 'AAAAACACATG': 'VVVVVVVVVVVVVVVC', 'AAAAAAAG': 'VVVVVVVVVVVVB?????==;', 'ACTGAG': 'VVVVVVVVVVL', 'AAGGCC': 'VVVVVVVVVVIG', 'AATAG': 'VVVVVVVVVN', 'AATAT': 'VVVVVVVVVIHE', 'AAAAAAAT': 'VVVVVVVVVVVVF@@@>>>>>>', 'CCCCCCG': 'VVVVVVVVVVVAAAA:', 'ACACACC': 'VVVVVVVVVVVVA', 'AAAAGCT': 'VVVVVVVVVVVHG', 'ACATATAT': 'VVVVVVVVVVVVGGCA@===<;;;;;7', 'AAAAGCC': 'VVVVVVVVVVVI', 'AAATAGAT': 'VVVVVVVVVVVVE', 'AAATATT': 'VVVVVVVVVVVMH', 'ACTGGC': 'VVVVVVVVVVJJ', 'AAAATATAT': 'VVVVVVVVVVVVVGGDA', 'AGCAGGG': 'VVVVVVVVVVVI', 'AACCAG': 'VVVVVVVVVVI', 'ACTGCGG': 'VVVVVVVVVVVVVV>', 'AAATATG': 'VVVVVVVVVVVLJ', 'AAAGGGG': 'VVVVVVVVVVVHGDA', 'AAGGCTAC': 'VVVVVVVVVVVVJ', 'ACACACAG': 'VVVVVVVVVVVVFFFA', 'AAAACAAACC': 'VVVVVVVVVVVVVVV?', 'ACACATATAT': 'VVVVVVVVVVVVVVA@===<<<<;;;;;;1', 'AAAAATATAC': 'VVVVVVVVVVVVVVB', 'AAGAAT': 'VVVVVVVVVVJIH', 'AAAAACCAT': 'VVVVVVVVVVVVVEE', 'AAGGC': 'VVVVVVVVVK', 'AAAAAGAAT': 'VVVVVVVVVVVVVVVVB', 'AGCCCCC': 'VVVVVVVVVVVII', 'AAGGG': 'VVVVVVVVVJJJJGGGG===;', 'ATCC': 'VVVVVVVKJJIIGG@>>><<888877766641111111000000,,,,,,,&', 'AAACACACAC': 'VVVVVVVVVVVVVVVV=', 'ACACACATAT': 'VVVVVVVVVVVVVVA>>>=', 'ACACCT': 'VVVVVVVVVVF', 'AAATC': 'VVVVVVVVVKKKH', 'AAGGT': 'VVVVVVVVVJ', 'AATTC': 'VVVVVVVVVML', 'AAAAATGTTTT': 'VVVVVVVVVVVVVVVC', 'AAAGT': 'VVVVVVVVVIIIH', 'ATATATC': 'VVVVVVVVVVVHH', 'AAACAGG': 'VVVVVVVVVVVH', 'AAAGACTT': 'VVVVVVVVVVVVH', 'AAAGC': 'VVVVVVVVVKK', 'AACTTCAGC': 'VVVVVVVVVVVVVI', 'AAAGG': 'VVVVVVVVVJIHHC', 'AAGGTGT': 'VVVVVVVVVVVI', 'AAATGGAT': 'VVVVVVVVVVVVVF', 'ACATGC': 'VVVVVVVVVVKKJ', 'ACC': 'VVVVVKJIHGDAA999885555444444222,', 'ACG': 'VVVVVNK', 'ACATAGGTAT': 'VVVVVVVVVVVVVVJ', 'AGCCTCC': 'VVVVVVVVVVVJD', 'AGGAGGC': 'VVVVVVVVVVVKF', 'ACAT': 'VVVVVVVMLLHFECA?>><999999999999999999**********', 'AGGAGGG': 'VVVVVVVVVVVIFB', 'ACT': 'VVVVVMLKII', 'AGAGG': 'VVVVVVVVVKKKKKD', 'AAGATT': 'VVVVVVVVVVK', 'AAAAAATGC': 'VVVVVVVVVVVVVEC', 'AGAGC': 'VVVVVVVVVL', 'AAAAAAT': 'VVVVVVVVVVVFDCCCCBA>>>>>=', 'AATGGT': 'VVVVVVVVVVIH', 'AGGGG': 'VVVVVVVVVIGGDDDDD>', 'CCCGG': 'VVVVVVVVVHC', 'AGGGC': 'VVVVVVVVVJJJJ', 'AAGATG': 'VVVVVVVVVVK', 'AATGGC': 'VVVVVVVVVVJ', 'AATGGG': 'VVVVVVVVVVI', 'ACAGGG': 'VVVVVVVVVVKJ', 'AAAAGATC': 'VVVVVVVVVVVVI', 'AGAT': 'VVVVVVVNLKFD@@???=:998766543322//......,,,,,,,,,,,,********((', 'ACCATCATC': 'VVVVVVVVVVVVVF====:7', 'AAATGATG': 'VVVVVVVVVVVVBBBBBBBB:', 'ACACCACC': 'VVVVVVVVVVVVVE', 'AAAATAAAT': 'VVVVVVVVVVVVVJCAAAAA>', 'AAAAGT': 'VVVVVVVVVVGG', 'AAAGAGG': 'VVVVVVVVVVVKI', 'AATATATAT': 'VVVVVVVVVVVVVVVV@', 'AAAAGG': 'VVVVVVVVVVHHHHE', 'AAAAACCCTTC': 'VVVVVVVVVVVVVVVAA', 'AAAAGC': 'VVVVVVVVVVKJF', 'AAAAAAATAAT': 'VVVVVVVVVVVVVVVBA', 'AAAATAAAC': 'VVVVVVVVVVVVVD', 'AAAAAGT': 'VVVVVVVVVVVGGG', 'AAACC': 'VVVVVVVVVLJJJJJJJ>', 'AAGTCC': 'VVVVVVVVVVJ', 'AAAAAGG': 'VVVVVVVVVVVHF', 'AAACT': 'VVVVVVVVVIIII', 'AAAAAGC': 'VVVVVVVVVVVIF', 'AACCTC': 'VVVVVVVVVVJ', 'AAAACTG': 'VVVVVVVVVVVJI', 'AAACTAT': 'VVVVVVVVVVVHHH', 'AAAAGGG': 'VVVVVVVVVVVHD', 'ACTGGCCATC': 'VVVVVVVVVVVVVVIE', 'AAAAATTT': 'VVVVVVVVVVVVIFFC', 'AAGAGGAGG': 'VVVVVVVVVVVVVFFB@', 'AATCAT': 'VVVVVVVVVVKE', 'AACAGC': 'VVVVVVVVVVK', 'AAATATAT': 'VVVVVVVVVVVVGGGE', 'AAAAATTC': 'VVVVVVVVVVVVGEE', 'AGATCAGG': 'VVVVVVVVVVVVVD', 'AC': "VVVOONLIFC>:86533220/.-,+++*****))))))))'''''''''''''''", 'AAAAAGTG': 'VVVVVVVVVVVVH', 'AG': 'VVVMMMLJHC?;:77442210--,+++*******))))))))))((&&&%', 'AGCAGG': 'VVVVVVVVVVLLG', 'AAAAAGTT': 'VVVVVVVVVVVVHF', 'ACCCTCC': 'VVVVVVVVVVVFD', 'AT': "VVVONNLIF@=:75310.++****)))''''''''''''''''", 'AAAATATT': 'VVVVVVVVVVVVG', 'AATAGG': 'VVVVVVVVVVK', 'AAATATATAT': 'VVVVVVVVVVVVVVVAAAA', 'AGAGAT': 'VVVVVVVVVVK', 'AAAATATG': 'VVVVVVVVVVVVG', 'AATAGT': 'VVVVVVVVVVLK', 'AACAATGAG': 'VVVVVVVVVVVVVML'}



def wrap(line, length=80):
    if len(line)<=length: return line
    if len(line)==0: return ""
    while length>=0 and line[length] != " ": length -= 1
    start = length
    while length>=0 and line[length] == " ": length -= 1
    if length>=0: return line[:length+1] + "\n" + wrap(line[start+1:])
    while length<len(line) and line[length] != " ": length += 1
    return line[:length] + "\n" + wrap(line[length+1:])


def help():
    print "Add homopolymer run length and tandem repeat length (HR and TR fields) to INFO column"
    print "Usage: vcf-add-hr.py [options]"
    print
    print " -h, --help         This help"
    print " -oF, --output F    Send output to file F (default stdout)"
    print " -iF, --input F     Read input from file F (default stdin)"
    print " -rF, --reference F Use reference fasta file; must be samtools faidx'ed (required)"
    print " -tN, --tandem N    Annotate for tandem repeat units of up to N bases (default 4)"
    print " --unstranded       Report unstranded homopolymer and tandem repeat unit (default: stranded)"
    print " --slippage         Annotate for potential slippage-associated indel"
    print " --indelrate        Annotate with model-based indel rate [experimental]"
    print " --hotspot          Annotate hotspots based on homopolymer/tandem repeat context"
    print " --repwindow N      Annotate hottest hotspot category in window of size N"
    print " --errorrate        Annotate with indel error rate, based on recent Illumina data"
    print " --palindrome N     Annotate with maximum palindromic match between ref and alt alleles"
    print " --transversion     Annotate transversions"
    print " --leftalign        Normalize (left-align) indel alleles"
    print " --region R         Annotate per-base in region R (e.g. chr1:1-1000)"
    print " --regionir N       Consider IH==2 | IR >= N a hotspot (default 20)"
    print " --gc N             Annotate for GC content in N bp windows"
    print " --nohr              Do NOT annotate homopolymers or tandems"
    print " -m, --mask         Annotate call by mask file (-r); implies --nohr"
    print " --mark F,L,D       Mark sites with info label L, description D; F lists sites as chrom-tab-pos pairs"
    print " -f C, --filter C   Filter call if mask file (-r) is not C; do not annotate homopolymers or tandems"
    print " -3                 Set default input version to 3.3"
    print " -xE, --ignore E    Ignore error E"
    print " -XE, --warn E      Treat E as warning"
    print
    print "Possible errors are:"
    v = vcf.VCF()
    print wrap( ", ".join([v.split(':')[0] for v in v._errors.values()]) )


def lcs(s, t, a=0, b=1e10):
    """ longest common substring, which includes at least one character from s[a:b].  Returns length, starting positions in s and t"""
    l0 = [0] * len(t)   # current row
    l1 = [0] * len(t)   # previous row
    z = 0               # lcs
    starts = -1
    startt = -1
    for i, sc in enumerate(s.upper()):
        for j, tc in enumerate(t.upper()):
            if sc == tc:
                if i==0 or j==0:
                    if i<b:
                        l0[j] = 1
                    else:
                        l0[j] = 0
                else:
                    if i<b or l1[j-1]>0:
                        l0[j] = l1[j-1] + 1
                    else:
                        l0[j] = 0
                if l0[j] >= z and i >= a:
                    if l0[j] > z or abs( startt + (z - len(t))//2 ) > abs( j-z+1 + (z - len(t)//2) ):
                        z = l0[j]
                        starts = i-z+1
                        startt = j-z+1
            else:
                l0[j] = 0
        l0, l1 = l1, l0
    #print "(lcs) matched:"
    #print s[:starts].lower() + s[starts:starts+z].upper() + s[starts+z:].lower()
    #print t[:startt].lower() + t[startt:startt+z].upper() + t[startt+z:].lower()
    return z, starts, startt
    

def get_max_palindrome(chrom, pos, fa, ref, alt, windowsize):
    """ Returns length of longest palindromic match between the ref and alt alleles, which overlaps
        the longest allele by at least 1 nucleotide.  This is symmetric for insertions vs deletions.
        Choose ref == alt == '' for a palindromic match of the reference sequence """
    
    seq = vcf.get_sequence(chrom, pos - windowsize, pos + windowsize + max(len(ref),len(alt)), fa)
    assert seq[ windowsize : windowsize + len(ref) ] == ref
    seq2 = seq[ :windowsize ] + alt + seq[ windowsize+len(ref) : ]
    if len(alt)>len(ref):
        # insertion
        #print "Comparing alt to -ref (ins); alt=",seq2[:windowsize+1].lower() + seq2[windowsize+1:windowsize+len(alt)].upper() + seq2[windowsize+len(alt):].lower()
        lng,strt1,strt2 = lcs( seq2, revcmp(seq), windowsize+1, windowsize + len(alt) )
        if strt2 > -1:
            return lng, pos - windowsize + (len(seq) - strt2 - lng)
        else:
            return lng, -1
    else:
        # deletion
        #print "Comparing ref to -alt (del); ref=",seq[:windowsize+1].lower() + seq[windowsize+1:windowsize+len(ref)].upper() + seq[windowsize+len(ref):].lower()
        lng,strt1,strt2 = lcs( seq, revcmp(seq2), windowsize+1, windowsize + len(ref) ) 
        if strt1 > -1:
            return lng, pos - windowsize + strt1
        else:
            return lng, -1


def microsat( seq, start, replen, direction ):
    minpos = maxpos = start
    while minpos>=0 and minpos+replen*direction>=0 and seq[minpos] == seq[minpos + replen*direction] and seq[minpos] != 'N': 
        minpos -= 1
    while maxpos<len(seq) and maxpos+replen*direction<len(seq) and seq[maxpos] == seq[maxpos + replen*direction] and seq[maxpos] != 'N': maxpos += 1
    length = maxpos - minpos + replen - 1
    if length < min(replen + min_matching_stretch, 2*replen):
        return 1,seq[start]
    return length, seq[minpos+1:minpos+replen+1]


def microsattelite( seq, start, tandem=4, backward=True ):
    maxlen = (-1,"")
    if backward: dirs = ((0,-1),(0,1),(-1,-1),(-1,1))
    else:        dirs = ((0,-1),(0,1))
    for m in range(1,tandem+1):
        for dx, direction in dirs:
            length, unit = microsat(seq, start+dx, m, direction)
            if length > maxlen[0]: maxlen = length,unit
    return maxlen


def get_hotspot_category( homopolsize, size, unit ):
    if ( homopolsize >= tandem_thresholds_high[1] or 
         (len(unit) < len(tandem_thresholds_high) and size >= tandem_thresholds_high[len(unit)])):
        return 2
    elif (homopolsize >= tandem_thresholds[1] or
          size >= tandem_thresholds[min(len(unit),len(tandem_thresholds)-1)]):
        return 1
    return 0


def get_repetitive_window(chrom, pos, fa, windowsize, tandem):
    """ returns hotspot category: 0 complex, 1 warm, 2 hotspot """

    seqwindow = max(20, tandem+10)
    reppos = pos - windowsize
    maxcategory = 0
    seq = vcf.get_sequence(chrom, pos - windowsize - seqwindow, pos + windowsize + seqwindow, fa)
    while reppos < pos + windowsize:
        t, repunitT = microsattelite( seq[reppos - seqwindow - (pos-windowsize-seqwindow) :
                                          reppos + seqwindow - (pos-windowsize-seqwindow)], 
                                      seqwindow, tandem=tandem, backward=True )
        
        cat = get_hotspot_category( 0, t, repunitT )

        if cat > maxcategory:
            maxcategory = cat
            if maxcategory == 2:
                break
        
        if t > seqwindow*0.33:
            seqwindow *= 2
            seq = vcf.get_sequence(chrom, pos - windowsize - seqwindow, pos + windowsize + seqwindow, fa)
        else:
            reppos += 1
    return maxcategory


def get_homopolymer_and_tandem(chrom, pos, fa, tandem, stranded, backward=True):
    window = max(20, tandem+10)
    while True:
        seq = vcf.get_sequence(chrom, pos-window, pos+window, fa)
        h, repunit1 = microsattelite( seq, window, tandem=1, backward=backward )
        t, repunitT = microsattelite( seq, window, tandem=tandem, backward=backward )
        if max(h,t) < window*0.33: 
            return h,t,normalize_repunit(repunit1, stranded),normalize_repunit(repunitT, stranded),seq
        window *= 2


def get_gcregion(chrom, pos, fa):
    window = 20
    while True:
        seq = vcf.get_sequence(chrom, pos-window, pos+window, fa).upper()
        idx1, idx2 = window, window
        while idx1 >= 0 and seq[idx1] in "GC": idx1 -= 1
        while idx2 < len(seq) and seq[idx2] in "GC": idx2 += 1
        if idx1>0 and idx2 < len(seq):
            size = max(0, idx2-idx1-1)
            return size, size, "G", "G", seq
        window *= 2



def get_indel_error_rate(tandemlength, tandemunit):
    if tandemunit.upper() in indel_error_model:
        model = indel_error_model.get(tandemunit.upper())
    elif len(tandemunit) in indel_error_model:
        model = indel_error_model.get(len(tandemunit))
    else:
        model = indel_error_model[1]
    rate = model[ min( len(model)-1, tandemlength-1 ) ]
    # convert to Phred
    return ord(rate) - 33


def get_indelrate(chrom, pos, fa):
    max_gap = 24
    max_match = 17
    window = 2*(max_gap+max_match+1)
    seq = vcf.get_sequence(chrom, pos-window, pos+window, fa)
    start = window
    rate = 1
    for gap in range(1,max_gap+1):
        rates = []
        for direction in (1,-1):
            # find longest match in indicated direction
            mmin = mmax = 0
            try:
                while mmax < max_match and seq[start + mmax] == seq[start + mmax + direction*gap]:
                    mmax += 1
            except:
                mmax -= 1
            try:
                while mmax-mmin-1 < max_match and seq[start + mmin] == seq[start + mmin + direction*gap]:
                    mmin -= 1
            except:
                mmin += 1
            matchsize = mmax-mmin-1
            thisrate = 0
            for m in range(1,min(matchsize+1,len(_r))):
                thisrate += _r[m] * _mu[gap]
            rates.append(thisrate)
        rate += max(rates)
    return int(0.5 + 10*math.log(rate)/math.log(10))
            
            

def get_gc(chrom, pos, fa, window):
    seq = vcf.get_sequence(chrom, pos-window//2, pos-window//2 + window, fa).upper()
    gc = seq.count('G') + seq.count('C')
    at = seq.count('A') + seq.count('T')
    if gc+at == 0: return 0
    return ( 1000*gc + (gc+at)//2 ) // (gc+at)


def getslippage(seq, chrom, pos, fa, ref, alt):
    seq = seq.upper()
    ref = ref.upper()
    alt = alt.upper()
    size = abs(len(ref)-len(alt))
    window = len(seq) // 2
    if (size+1)*2 >= window:
        window = 2*(size+2)
        seq = vcf.get_sequence(chrom, pos-window, pos+window, fa)
    pos = window-1  # compensate for 1-bp lead for VCF indel calls
    # if the variant is very long, or in case of a reference mismatch, return False, 0
    if max(len(ref),len(alt)) > window or seq[pos:pos+len(ref)] != ref:
        #print "** Long, or ref mismatch:",seq[pos:pos+len(ref)], ref
        return False, 0 
    # if the variant is complex, return False, 0
    if ref[:min(len(ref),len(alt))] != alt[:min(len(ref),len(alt))]:
        #print "** Complex"
        return False, 0
    # build insertion allele
    if len(alt) < len(ref):
        insallele = seq
        delallele = seq[:pos] + alt + seq[pos + len(ref):]
    else:
        insallele = seq[:pos] + alt + seq[pos + len(ref):]
        delallele = seq
    # check for slippage
    slippage = False
    for dx in range(-size,0) + range(1,size+1):
        for i in range(size):
            # account for 1-bp lead
            if insallele[1+pos+dx+i] != insallele[1+pos+i]:
                break
        else:
            slippage = True
            break
    # calculate length of direct repeat of the inserted sequence fragment
    dr = []
    left,right = 1+pos, pos+size
    while left+size < len(insallele) and left <= pos+size and insallele[left] == insallele[left+size]:
        left += 1
    while right-size >= 0 and right >= left and insallele[right-size] == insallele[right]:
        right -= 1
    direct_repeat = left - right + size - 1    
    return slippage, direct_repeat


def revcmp( unit ):
    return ''.join(reversed([{'A':'T','T':'A','C':'G','G':'C'}.get(c,'N') for c in unit.upper()]))


def normalize_repunit( unit, stranded ):
    if not stranded:
        unit2 = revcmp( unit )
    else:
        unit2 = unit
    return sorted([ unit[i:]+unit[:i] for i in range(len(unit)) ] + [ unit2[i:]+unit2[:i] for i in range(len(unit)) ])[0]


def readmark( filename ):
    data = {}
    for line in filez.open(filename,'r'):
        chrom, pos = line[:-1].split('\t')
        if chrom not in data: data[chrom] = []
        data[chrom].append(int(pos) - 1)   # convert into 0-based coords
    for chrom in data.keys():
        data[chrom] = set(data[chrom])
    return data


class Localcalls:

    def __init__(self, v, infile, windows):
        self.v = vcf.VCF( v )
        self.vcfstream = self.v.parse( infile, parseGenotypes=False )
        self.windows = windows
        self.maxwindow = max(windows)
        self.buf = []
        self.chromdict = {}

    def move(self, chrom, pos):
        
        # read records until we have enough
        while (len(self.buf)==0 or 
               self.chromdict[self.buf[-1][0]] < self.chromdict.get(chrom,1e10) or
               (self.chromdict[self.buf[-1][0]] == self.chromdict[chrom] and self.buf[-1][1] < pos + self.maxwindow ) ):

            try:
                data = self.vcfstream.next()
            except StopIteration:
                break

            if data['chrom'] not in self.chromdict:
                self.chromdict[data['chrom']] = len(self.chromdict)

            self.buf.append( (data['chrom'], data['pos']) )

        # remove records we don't need
        while len(self.buf)>0 and (self.chromdict[self.buf[0][0]] < self.chromdict[chrom] or self.buf[0][1] < pos - self.maxwindow):
            del self.buf[0]

        # get position in current chromosome
        positions = [ rec[1] for rec in self.buf if rec[0] == chrom ]

        return [ bisect.bisect(positions, pos+w) - bisect.bisect(positions, pos-w)
                 for w in self.windows ]
    

def main():
    infile = sys.stdin
    infile2 = None
    outfile = sys.stdout
    inversion = 40
    outversion = 40
    tandem = 4
    reference = None
    v = vcf.VCF( _fastGT = False )
    context = 0
    stranded = True
    mask = False
    slippage = False
    indelrate = False
    hotspot = False
    errorrate = False
    region = False
    regionir = 20
    transversion = False
    localcalls = None
    gc = 0
    repwindowsize = 0
    palindrome = 0
    palpos = None
    nohr = False
    mark = None
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ho:i:x:X:r:3mf:t:F", ["help","output=","input=","ignore=","warn=","reference=","tandem=","addcontext=","destranded","mask","filter=","slippage","gc=","indelrate","hotspot","leftalign","outversion=","region=","errorrate","repwindow=","palindrome=","transversion","localcalls=","nohr","mark=","regionir="])
    except:
        help()
        raise
    for o, a in opts:
        if o in ["-h","--help"]:
            help()
            sys.exit()
        elif o in ["-o","--output"]:
            outfile = open(a,'w')
        elif o in ["-i","--input"]:
            infile = filez.open(a,'r')
            infile2 = filez.open(a,'r')  # for localcalls
        elif o == "-3":
            inversion = 33
        elif o == "--outversion":
            outversion = int(a)
        elif o in ["-x","--ignore"]:
            v.ignoreerror(a)
        elif o in ["-X","--warn"]:
            v.warnerror(a)
        elif o in ["-r","--reference"]:
            reference = a
        elif o in ["-t","--tandem"]:
            tandem = int(a)
        elif o in ["-m","--mask"]:
            mask = "mask"
        elif o in ["-f","--filter"]:
            mask = a  
        elif o in ["--leftalign"]:
            v._leftalign = True
        elif o in ["-F"]:
            v._fastGT = True
        elif o in ["--localcalls"]:
            localcalls = map(int, a.split(','))
        elif o in ["--addcontext"]:
            context = int(a)
        elif o in ["--destranded"]:
            stranded = False
        elif o in ["--slippage"]:
            slippage = True
        elif o in ["--errorrate"]:
            errorrate = True
        elif o in ["--gc"]:
            gc = int(a)
        elif o in ["--indelrate"]:
            indelrate = True
        elif o in ["--palindrome"]:
            palindrome = int(a)
        elif o in ["--repwindow"]:
            repwindowsize = int(a)
        elif o in ["--hotspot"]:
            hotspot = True
        elif o in ["--transversion"]:
            transversion = True
        elif o in ["--nohr"]:
            nohr = True
        elif o in ["--mark"]:
            mark = a.split(',')
            assert len(mark) == 3
        elif o in ["--region"]:
            region = True
            chrom,startend = a.split(':')
            start, end = map(int, startend.split('-'))
        elif o in ["--regionir"]:
            regionir = int(a)
            
    if not reference:
        print "Reference required"
        help()
        sys.exit()

    if localcalls and not infile2:
        print "Cannot use stdin with --localcalls"

    # open reference
    fa = pysam.Fastafile( reference )
    if not mask: v.setreference(fa)
    
    # annotate region
    if region:
        if do_palindromes:
            collect_palindromes(fa, chrom, start, end, tandem, palindrome)
        else:
            annotate_region(fa, chrom, start, end, tandem, indelrate, regionir)
        return

    # read mark data
    if mark: markdata = readmark( mark[0] )

    # process data
    v.setversion(inversion)
    vcfstream = v.parse( infile , parseGenotypes=True)

    if localcalls:
        theLocalcalls = Localcalls( v, infile2, localcalls )

    # instantiate vcfout from v, to include header and definitions.
    vcfout = vcf.VCF(v)
    vcfout.setversion(outversion)
    vcfout.getheader().append( ("source",' '.join(["vcf-add-hr.py"]+sys.argv[1:])) )
    if stranded: destranded = " (stranded)"
    else:        destranded = " (unstranded)"
    if mask:
        if mask == "mask": vcfout.getinfo()['MASK'] = vcf.FORMAT('MASK',vcfout.NT_NUMBER,1,"Character","Mask",".")
        else:              vcfout.getfilter()['mask'] = vcf.FORMAT('mask',vcfout.NT_NUMBER,0,"Flag","Position masked",".")
    else:
        if not nohr:
            vcfout.getinfo()['HR'] = vcf.FORMAT('HR',vcfout.NT_NUMBER,1,"Integer","Homopolymer run length",-1)
            vcfout.getinfo()['HU'] = vcf.FORMAT('HU',vcfout.NT_NUMBER,1,"String","Homopolymer run unit%s" % destranded,-1)
            vcfout.getinfo()['TR'] = vcf.FORMAT('TR',vcfout.NT_NUMBER,1,"Integer","Tandem repeat run length (bp)",-1)
            vcfout.getinfo()['TU'] = vcf.FORMAT('TU',vcfout.NT_NUMBER,1,"String","tandem repeat run unit%s" % destranded,-1)
        if mark:
            vcfout.getinfo()[ mark[1] ] = vcf.FORMAT( mark[1], vcfout.NT_NUMBER,0,"Flag", mark[2], -1)
        if slippage:
            vcfout.getinfo()['SL'] = vcf.FORMAT('SL',vcfout.NT_NR_ALLELES,0,"Character","Indel appears to have been caused by a polymerase slippage event",".")
            vcfout.getinfo()['DR'] = vcf.FORMAT('DR',vcfout.NT_NR_ALLELES,0,"Integer","Length of direct repeat copy of long allele",-1)
        if palindrome > 0:
            vcfout.getinfo()['PAL'] = vcf.FORMAT('PAL',vcfout.NT_NUMBER,1,"Integer","Length of palindromic match between REF and first ALT allele",-1)
        if palindrome < 0:
            vcfout.getinfo()['PAL'] = vcf.FORMAT('PAL',vcfout.NT_NUMBER,1,"Integer","Length of maximum palindromic match on reference strand",-1)
        if hotspot:
            vcfout.getinfo()['IH'] = vcf.FORMAT('IH',vcfout.NT_NUMBER,1,"Character","Indel hotspot",".")
        if transversion:
            vcfout.getinfo()['TV'] = vcf.FORMAT('TV',vcfout.NT_NUMBER,1,"Integer","1 for transversions; 0 for transitions",".")
        if localcalls:
            vcfout.getinfo()['LC'] = vcf.FORMAT('LC',vcfout.NT_NUMBER,len(localcalls),"Integer","Number of local calls, in window(s) of size(s) %s" % ','.join(map(str,localcalls)), ".")
        if repwindowsize>0:
            vcfout.getinfo()['IHW'] = vcf.FORMAT('IHW',vcfout.NT_NUMBER,1,"Character","Indel hotspot in size-%s nt window" % repwindowsize, ".")
        if indelrate:
            vcfout.getinfo()["IR"] = vcf.FORMAT('IR',vcfout.NT_NUMBER,1,"Integer","Model-based location-specific indel rate, expressed as phred score indicating relative increase above base rate; e.g. 10=10x, 20=100x increase",-1)
        if errorrate:
            vcfout.getinfo()["IER"] = vcf.FORMAT('IER',vcfout.NT_NUMBER,1,"Integer","Estimated indel error rate, expressed as a Phred score of errors per read per repeat locus",-1)
        if gc>0:
            gclabel = 'GC'+str(gc)
            vcfout.getinfo()[gclabel] = vcf.FORMAT(gclabel,vcfout.NT_NUMBER,1,"Integer","GC content fraction in %s bp windows, times 1000" % gc,".")
    vcfout.writeheader( outfile )

    for data in vcfstream:
        chrom, pos = data['chrom'], data['pos']
        snp = True
        if mark:
            if pos in markdata.get(chrom, []):
                data['info'][ mark[1] ] = []
        # for indels, skip leading base
        if len(data['ref']) != 1 or sum(len(a) for a in data['alt']) != len(data['alt']):
            pos += 1
            snp = False
        if mask:
            m = fa.fetch(data['chrom'],pos,pos+1)
            if mask == "mask": data['info']['MASK'] = [m]
            elif m != mask: data['filter'].append('mask')
        else:
            if (not nohr) or errorrate or hotspot:
                homopolymer, tandemlen, homopolymerunit, tandemunit, seq = get_homopolymer_and_tandem(chrom, pos, fa, tandem, stranded)
            if not nohr:
                data['info']['HR'] = [homopolymer]
                data['info']['HU'] = [homopolymerunit]
                data['info']['TR'] = [tandemlen]
                data['info']['TU'] = [tandemunit]
            if errorrate:
                rate = get_indel_error_rate( tandemlen, tandemunit )
                data['info']['IER'] = [ rate ]
            if slippage:
                slip = [ getslippage(seq, chrom, pos, fa, data['ref'], alt) for alt in data['alt'] ]
                data['info']['SL'] = [ "NY"[sl[0]] for sl in slip ]
                data['info']['DR'] = [ sl[1] for sl in slip ]
            if repwindowsize>0:
                category = get_repetitive_window(chrom, pos, fa, repwindowsize, tandem)
                data['info']['IHW'] = [str(category)]
            if palindrome!=0:
                # negative window sizes signals the use of allele-independent palindromes (uses reference only)
                if palindrome > 0:
                    pallen, palpos = get_max_palindrome(chrom, data['pos'], fa, data['ref'], data['alt'][0], palindrome)
                else:
                    pallen, palpos = get_max_palindrome(chrom, data['pos'], fa, "", "", -palindrome)
                data['info']['PAL'] = [pallen]
            if transversion:
                if snp:
                    if data['ref']+data['alt'][0] in ["AG","GA","TC","CT"]:
                        data['info']['TV'] = [0]
                    else:
                        data['info']['TV'] = [1]
            if context!=0: 
                seq = fa.fetch(data['chrom'],data['pos']-abs(context),data['pos']+abs(context))
                # upper-case the palindrome match positions on the reference
                if palpos:
                    seq = seq.lower()
                    if palpos > -1:
                        posl, posr = max(0,palpos - (data['pos'] - abs(context))), min(palpos + pallen - (data['pos'] - abs(context)), len(seq))
                        seq = seq[:posl] + seq[posl:posr].upper() + seq[posr:]
                if context>0:
                    data['id'] = seq[:context] + "." + seq[context:] + ":" + data['id']
                else:
                    data['id'] = seq + ":" + data['id']
            if gc>0:
                ggcc = get_gc(chrom, pos, fa, gc)
                data['info'][gclabel] = [ggcc]
            if localcalls:
                data['info']['LC'] = theLocalcalls.move( chrom, data['pos'])
            if indelrate:
                indelrateest = get_indelrate(chrom, pos, fa)
                data['info']['IR'] = [indelrateest]
            if hotspot:
                cat = get_hotspot_category( homopolymer, tandemlen, tandemunit )
                data['info']['IH'] = [str(cat)]

        vcfout.write_data( outfile, data )


def collect_palindromes(fa, chrom, start, end, tandem, palindrome):

    palhist = {}
    pos = start
    rep_skip = 2
    skip_until = -1
    while pos < end:

         homopolymer, tandemlen, homopolymerunit, tandemunit, seq = get_homopolymer_and_tandem(chrom, pos + rep_skip, fa, tandem, 
                                                                                               False, backward=False)
         if seq.count('N') == len(seq):
             pos += len(seq)
             continue

         if (homopolymer >= tandem_thresholds[1] or
             (len(tandemunit)<len(tandem_thresholds) and tandemlen >= tandem_thresholds[len(tandemunit)]) or
             len(tandemunit)>=len(tandem_thresholds)):
             
             skip_until = max( skip_until, pos + max(homopolymer, tandemlen) + 2*rep_skip)

         if pos < skip_until:
             pos += 1
             continue
        
         pallen, palpos = get_max_palindrome(chrom, pos, fa, "", "", palindrome)

         palhist[pallen] = palhist.get(pallen,0) + 1

         pos += 1

    # done
    for key in sorted(palhist.keys()):
        print "%s\t%s" % (key, palhist[key])


def annotate_region(fa, chrom, start, end, tandem, indelrate, regionir):
    """ Annotate regions by homopolymer and tandem repeat content.  Concatenate inexact tandem repeats consisting of the same unit """

    global curseq
    global notmasked
    curseq = ">%s:%s-%s" % (chrom,start,end)

    notmasked = 0

    pos = start-1
    oldoutput = ["", -1, -1, "", -1, ""]

    while pos < end:

        if do_get_gcregion:
            homopolymer, tandemlen, homopolymerunit, tandemunit, seq = \
                get_gcregion(chrom, pos, fa)
        else:
            homopolymer, tandemlen, homopolymerunit, tandemunit, seq = \
                get_homopolymer_and_tandem(chrom, pos, fa, tandem, False, backward=False)

        if seq.count('N') == len(seq):
            if do_mask:
                output_fasta(seq[ len(seq)/2 : ])
            pos += len(seq)/2
            continue

        output = [chrom, pos, homopolymer, homopolymerunit, tandemlen, tandemunit]

        if indelrate or do_mask:
            indelrateest = get_indelrate(chrom, pos, fa)
            output += [indelrateest]

        # concatenate or output
        if do_mask:
            if do_get_gcregion:
                cat = (homop > 6)
                regionir = 0
            else:
                cat = get_hotspot_category( homopolymer, tandemlen, tandemunit )
            if cat == 2 or indelrateest >= regionir:
                size = max(homopolymer, tandemlen)
                output_fasta( seq[ len(seq)/2: len(seq)/2 + size ].lower() )
                pos += max(size,1)
            else:
                output_fasta( seq[ len(seq)/2 ].upper() )
                pos += 1
        else:
            if homopolymer >= 6 or tandemlen >= 7 or (get_gcregion and homopolymer >= 6):
                if oldoutput[1] + oldoutput[4] == output[1] and oldoutput[5] == output[5]:
                    oldoutput[4] += output[4]
                else:
                    if oldoutput[1]>0: print (("%s\t" * len(oldoutput)) % tuple(oldoutput))[:-1]
                    oldoutput = output
                pos += max(output[2],output[4])
            else:
                pos += 1

    # done
    if do_mask:
        output_fasta()
        output_fasta()
    elif oldoutput[1]>0:
        print (("%s\t" * len(oldoutput)) % tuple(oldoutput))[:-1]


def output_fasta( seq=None ):

    global curseq
    global notmasked

    if curseq.startswith(">"):
        print curseq
        curseq = ""

    if seq == None:
        seq = ""
    if len(seq) > 1 and notmasked > 0 and notmasked <= mask_merge:
        curseq = curseq[:-notmasked] + curseq[-notmasked:].lower()
    curseq += seq
    if len(seq) != 1:
        notmasked = 0
    else:
        notmasked += 1
    if len(curseq) >= 80 + mask_merge:
        print curseq[:80]
        curseq = curseq[80:]
    elif len(seq) == 0:
        if len(curseq) > 0:
            toprint = max( len(curseq), 80 )
            print curseq[:toprint]
            curseq = curseq[toprint:]

if __name__ == "__main__": main()



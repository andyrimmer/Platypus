"""
module containing various classes and functions for use in generating
and processing haplotypes.
"""

from __future__ import division

import cython
import logging
import math
import random

cimport cython
cimport variant
cimport calign
cimport fastafile
cimport samtoolsWrapper
cimport cerrormodel

from calign cimport hash_sequence, hash_sequence_multihit, hashReadForMapping
from calign cimport mapAndAlignReadToHaplotype
from fastafile cimport FastaFile
from samtoolsWrapper cimport AlignedRead
from samtoolsWrapper cimport cAlignedRead
from samtoolsWrapper cimport Read_IsReverse
from samtoolsWrapper cimport Read_IsPaired
from samtoolsWrapper cimport Read_IsProperPair
from samtoolsWrapper cimport Read_IsDuplicate
from samtoolsWrapper cimport Read_IsUnmapped
from samtoolsWrapper cimport Read_MateIsUnmapped
from samtoolsWrapper cimport Read_MateIsUnmapped
from samtoolsWrapper cimport Read_IsQCFail
from samtoolsWrapper cimport Read_IsReadOne
from samtoolsWrapper cimport Read_IsSecondaryAlignment
from variant cimport Variant
from calign cimport hash_nucs,hash_size

###################################################################################################

logger = logging.getLogger("Log")

###################################################################################################

cdef double PI = math.pi
cdef double mLTOT = -0.23025850929940459    # Minus log ten over ten

cdef extern from "math.h":
    double exp(double)
    double log(double)
    double log10(double)
    double fabs(double)
    double sqrt(double)
    double pow(double, double)

cdef extern from "stdlib.h":
    void free(void *)
    void *malloc(size_t)
    void *calloc(size_t,size_t)
    void *realloc(void *,size_t)
    void *memset(void *buffer, int ch, size_t count )

###################################################################################################

# Some nasty global variables
cdef list per_base_indel_errors = [2.9e-5, 2.9e-5, 2.9e-5, 2.9e-5, 4.3e-5, 1.1e-4, 2.4e-4, 5.7e-4, 1.0e-3, 1.4e-3] + [ 1.4e-3 + 4.3e-4*(n-10) for n in range(11,50) ]

# homopolymer indel error model
cdef bytes homopolq = bytes(''.join([chr(int(33.5 + 10*log( (idx+1)*q )/log(0.1) )) for idx,q in enumerate(per_base_indel_errors)]))

# same model, represented in a more general way.  Commented out for now
#cdef dict indel_error_model = {'nonrepetitive':'X',
#                               1:'NKJHFA=854210/.-,,+**))(((\'\'\'&&&%%%$$$$#####"""""'}

# model obtained from 10M reads from AW_SC_4654.bam, using an early version of indelerrormodel.py (updated the default extrapolation exponents)

default_indel_error_model = {'AAGG': 'WWWWWWWKJHFEC<<;98665320/-,*)(&%#"!', 1: 'WWWKIFB<83/+(%"!!!!!!!!!!!!!!', 2: "WWWKKKJJFB=:86531/.,*('%#!!!!!!!!!!!!!!!!!!!!!!", 3: 'WWWWWLLKIGEDB@====;;53220.,+)\'%#"!!!', 4: "WWWWWWWKJHHFFBBBB>>=<<;764310.-+*)'&$#!!!!!!!!!!!!!!!!!!!!!!!!", 5: 'WWWWWWWWWIIIFEA@???><<<;9865320/-,+)(&%#"!', 6: 'WWWWWWWWWWIIIIFDCA@>=;:87754210.-+*(\'%$"!', 7: 'WWWWWWWWWWWIHHEEA?>>=;:875421/.,+*(\'%$"!', 8: 'WWWWWWWWWWWWIFEECB@?;99865320/-,*)\'&%#"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', 9: 'WWWWWWWWWWWWWHFFFDB=;:875421/.-+*(\'%$"!', 10: "WWWWWWWWWWWWWWHFFCB??><;99764310.-+*)'&$#!!", 11: "WWWWWWWWWWWWWWWGGEEBB?=:9764310.-,*)'&$#!!", 'AGC': 'WWWWWMMKKKGFDB@?=;986420/-+)(&$"!', 'AAG': 'WWWWWMMLJHGECA@><:975320.,+)\'%#"!', 'ACCC': 'WWWWWWWHHHGEDBA?><;:875421/.,+)(\'%$"!', 'AAAAT': 'WWWWWWWWWIHFEDBA?><;9865321/.,+)(&%#"!', 'AAT': 'WWWWWKJIHGEC@><;;9775420.-+)\'%$"!!!!!!!!!!!', 'AAAC': "WWWWWWWJJHGFDCA@>==;:888754310.-+*('%$#!!", 'AAGT': "WWWWWWWLGFDCA@>=;:9764310.-+*('&$#!!", 'AAAG': 'WWWWWWWJJFEEA@?=<:9765320/-,*)\'&$#"!', 'ACT': 'WWWWWMJJIGEDB@><;975420.-+)\'&$"!', 'ACCT': "WWWWWWWJIGFDCA@?=<:9764310.-,*)'&$#!!", 'CCG': 'WWWWWIGECA@><:975320.,+)\'%$"!', 'AAAT': 'WWWWWWWJIHFEECCBA@@>=;;:875421/.,+*(\'%$"!', 'AAGC': 'WWWWWWWKFEDBA?><;9865421/.,+)(&%#"!', 'AAAAAATAC': "WWWWWWWWWWWWWA>>=;:9764310.-+*('&$#!!", 'AAAAATAT': "WWWWWWWWWWWWVTA@>=<:9764310.-+*)'&$#!!", 'AGGC': 'WWWWWWWJIHFEDBA?><;9865321/.,+)(&%#"!', 'AACC': "WWWWWWWKJIGFECB@?=<:9764320/-,*)'&$#!!", 'AGCC': "WWWWWWWJIHDCA@>=;:9764310.-+*('&$#!!", 'AAAAAAATAC': 'WWWWWWWWWWWWWW>>><;9865321/.,+)(&%#"!', 'AAAAAAT': 'WWWWWWWWWWWD@@?=<:9765320/-,*)\'&$#"!', 'AAAAATAC': "WWWWWWWWWWWWBA@?=;:9764310.-+*('&$#!!", 'ACG': 'WWWWWKIGEDB@>=;975420.-+)\'&$"!', 'ACCTCC': 'WWWWWWWWWWVHHFECBA?><;9865320/.,+)(&%#"!', 'AACT': 'WWWWWWWHGEDBA?><;9865421/.,+)(&%#"!', 'ACAG': 'WWWWWWWMKJHGEDCA@>=;:8754210.-+*(\'%$"!', 'AATG': "WWWWWWWMKJIGFDCA@>=;:8764310.-+*('%$#!!", 'AAAAAGAAAAG': 'WWWWWWWWWWWWWWWVT<;9865321/.,+)(&%#"!', 'AC': 'WWWNMLIGGC=987631/-,*(&$#!!!!!!!!!!!!!!!!!!!!!!', 'AAAAAG': 'WWWWWWWWWWFBA?><;9865321/.,+)(&%#"!', 'AG': 'WWWIIIIIDD?=9774420-,*(&%#!!!!!!!!!', 'AGGG': 'WWWWWWWHHEECB@?=<:9765320/-,*)\'&%#"!', 'AAACAGAC': 'WWWWWWWWWWWWVA?><;9865421/.,+)(&%$"!', 'CG': 'WWWMJHGECA@><:975320.,+)\'%#"!', 'AAAAAAG': 'WWWWWWWWWWWBA?><;9865320/.,+)(&%#"!', 'CCCG': "WWWWWWWHGFDCB@?=<:9764310/-,*)'&$#!!", 'AT': "WWWMMMKKE?<:7420/-+('%#!!", 'ACCCC': "WWWWWWWWWECB@?=<:9764320/-,*)'&$#!!", 'ACATATATAT': "WWWWWWWWWWWWWWV>=;:8764310.-+*('%$#!!", 'ACC': "WWWWWKKIHEDCA?=<:86531/.,*('%#!!", 'AATAT': 'WWWWWWWWWVCA@>=;:8754310.-+*(\'%$"!!', 'AGCCTC': "WWWWWWWWWWJIGFECB@?=<:9764320/-,*)'&$#!!", 'AGGCGG': "WWWWWWWWWWIIGFDCA@>=;:8764310.-+*('%$#!!", 'AAAAAT': 'WWWWWWWWWWHFEDBA?><;9865321/.,+)(&%#"!', 'AATT': 'WWWWWWWKIIHFECB@?=<;9865320/-,*)(&%#"!', 'C': "WWWJHD<74-*'$!", 'ACACACATAT': 'WWWWWWWWWWWWWW;:8754210.-+*(\'%$"!', 'AGGCCG': 'WWWWWWWWWWHGFDCA@>=;:8754310.-+*(\'%$"!', 'AAAATAC': "WWWWWWWWWWW@?=<:9764320/-,*)'&$#!!", 'AAAAGAAAG': 'WWWWWWWWWWWWWVT?><;9865320/.,+)(&%#"!', 'AAAAAC': 'WWWWWWWWWWEDBA?>=;:875421/.,+*(\'%$"!', 'ACAT': 'WWWWWWWLKKIHFECB@?=<:9865320/-,*)\'&%#"!', 'ATC': 'WWWWWMMKHFDCA?=<:86431/-,*(&%#!!', 'ATCCC': 'WWWWWWWWWECB@?=<:9765320/-,*)\'&$#"!', 'AGG': 'WWWWWKKJHFDB@?=;986421/-+*(&$#!', 'AAAAAAAATAC': 'WWWWWWWWWWWWWWW9765320/-,*)\'&$#"!', 'AAAAAAAG': 'WWWWWWWWWWWWV=<:9765320/-,*)\'&$#"!', 'AAAAC': 'WWWWWWWWWFECCB@?==<:::9765320/-,*)\'&$#"!', 'ACTG': 'WWWWWWWMKJHGEDBA?>=;:875421/.,+*(\'%$"!', 'AAAAG': 'WWWWWWWWWGFDCB>;:875421/.,+*(\'%$"!', 'ATCC': 'WWWWWWWLKIHFECB@?=<:9865320/-,*)\'&%#"!'}

# oth.errormodel.20111011.txt, in Research/1000GIndels/2011_09_27-error-rates/

#default_indel_error_model = { 1: 'DDDDDC?:620//.--,' ,
#                                2: 'DDDDDDDC?>:86442210///..-,,,,++++****)))))))))))(((%' ,
#                                3: 'DDDDDDDDDB>>=:::8876554331100///--------***)' ,
#                                4: 'DDDDDDDDDC===<888866665533222110//.------,,++++++*****)' ,
#                                5: 'DDDDDDDDDCBA=9999766666555333322........,,,,,,,,,)' ,
#                                6: "DDDDDDDDDDDB@?;:7777765555553322222211....---,,**(('''&" ,
#                                7: 'DDDDDDDDDDDDD?>=:8755555553333321' ,
#                                8: 'DDDDDDDDDDDDDDD==<:766444432221' ,
#                                9: 'DDDDDDDDDDDDDDDDD;;:7665422222222222110' ,
#                                10: 'DDDDDDDDDDDDDDDDDDD988764443222220/' ,
#                                11: 'DDDDDDDDDDDDDDDDDDDDD66644333000.' ,
#                                12: 'DDDDDDDDDDDDDDDDDDDDDDD6655554331111110..,' ,
#                                13: 'DDDDDDDDDDDDDDDDDDDDDDDDD55555432220' ,
#                                14: 'DDDDDDDDDDDDDDDDDDDDDDDDDDD6555432000000/' ,
#                                15: 'DDDDDDDDDDDDDDDDDDDDDDDDDDDDD66654' ,
#                                16: 'DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD7665220' ,
#                                17: 'DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD888842' ,
#                                18: 'DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD55555543' ,
#                                19: 'DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD7776' ,
#                                20: 'DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD887441111110/-+*)' ,
#                                21: 'DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD7742' ,
#                                22: 'DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD755443' ,
#                                23: 'DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD4' ,
#                                24: 'DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD652' ,
#                                'A': 'DDDDDC@:620//.--,' ,
#                                'AAAAAAAC': 'DDDDDDDDDDDDDDD5555554111111/' ,
#                                'AAAAAAAG': 'DDDDDDDDDDDDDDD5554222111110' ,
#                                'AAAAAAAT': 'DDDDDDDDDDDDDDD;;;;:876' ,
#                                'AAAAAAATT': 'DDDDDDDDDDDDDDDDD97' ,
#                                'AAAAAAC': 'DDDDDDDDDDDDD777776444444443321' ,
#                                'AAAAAACC': 'DDDDDDDDDDDDDDD766521' ,
#                                'AAAAAAG': 'DDDDDDDDDDDDD87743321111110' ,
#                                'AAAAAAGAAAAG': 'DDDDDDDDDDDDDDDDDDDDDDD2' ,
#                                'AAAAAAGAAAG': 'DDDDDDDDDDDDDDDDDDDDD320' ,
#                                'AAAAAAGG': 'DDDDDDDDDDDDDDD87521' ,
#                                'AAAAAAT': 'DDDDDDDDDDDDD>>>==;9999875' ,
#                                'AAAAAATAC': 'DDDDDDDDDDDDDDDDD88654' ,
#                                'AAAAAATT': 'DDDDDDDDDDDDDDD<<<:87' ,
#                                'AAAAAC': 'DDDDDDDDDDD<<<::77777777777763333310//-' ,
#                                'AAAAACC': 'DDDDDDDDDDDDD;;9865' ,
#                                'AAAAAG': 'DDDDDDDDDDD;:8643333211111110' ,
#                                'AAAAAGAAAAG': 'DDDDDDDDDDDDDDDDDDDDD430' ,
#                                'AAAAAGAAAG': 'DDDDDDDDDDDDDDDDDDD5522/' ,
#                                'AAAAAGAAAGAAAG': 'DDDDDDDDDDDDDDDDDDDDDDDDDDD4421' ,
#                                'AAAAAGAG': 'DDDDDDDDDDDDDDD;8' ,
#                                'AAAAAGG': 'DDDDDDDDDDDDD;;886' ,
#                                'AAAAAT': 'DDDDDDDDDDD???==;;;;;::752' ,
#                                'AAAAATAC': 'DDDDDDDDDDDDDDD<<7665' ,
#                                'AAAAATT': 'DDDDDDDDDDDDD=99997' ,
#                                'AAAAC': 'DDDDDDDDD?>>=:::::99999866444411/////...,,,,,,,+*&' ,
#                                'AAAACC': 'DDDDDDDDDDD@@>;:7' ,
#                                'AAAAG': 'DDDDDDDDD>><7444431111111110000.,' ,
#                                'AAAAGAAAG': 'DDDDDDDDDDDDDDDDD7652221' ,
#                                'AAAAGAAAGAAAG': 'DDDDDDDDDDDDDDDDDDDDDDDDD2' ,
#                                'AAAAGAAG': 'DDDDDDDDDDDDDDD<75' ,
#                                'AAAAGAG': 'DDDDDDDDDDDDD=:8' ,
#                                'AAAAGC': 'DDDDDDDDDDD@??:7' ,
#                                'AAAAGG': 'DDDDDDDDDDD??>:6' ,
#                                'AAAAT': 'DDDDDDDDDCCA@===<<<;;;885543331/.' ,
#                                'AAAATT': 'DDDDDDDDDDD@@@?:' ,
#                                'AAAC': 'DDDDDDDCBA<<<<::::8877664433211///..,' ,
#                                'AAACAC': 'DDDDDDDDDDDC==;64' ,
#                                'AAACACAC': 'DDDDDDDDDDDDDDD8853/.' ,
#                                'AAACACACAC': 'DDDDDDDDDDDDDDDDDDD540/-+' ,
#                                'AAAG': "DDDDDDDAA@88883333111111000000....----,,,,,+++***))))))((('" ,
#                                'AAAGAAAGAAAGAG': 'DDDDDDDDDDDDDDDDDDDDDDDDDDD4' ,
#                                'AAAGAAAGAG': 'DDDDDDDDDDDDDDDDDDD8731' ,
#                                'AAAGAAG': 'DDDDDDDDDDDDD==9876421' ,
#                                'AAAGAAGG': 'DDDDDDDDDDDDDDD=:870' ,
#                                'AAAGAG': 'DDDDDDDDDDDA==;:8887' ,
#                                'AAAGAGAG': 'DDDDDDDDDDDDDDD;9988654' ,
#                                'AAAGAGAGAG': 'DDDDDDDDDDDDDDDDDDD9876' ,
#                                'AAAGAGAGAGAGAG': 'DDDDDDDDDDDDDDDDDDDDDDDDDDD1/,' ,
#                                'AAAGG': 'DDDDDDDDDCCA=;;;;:87766555432' ,
#                                'AAAGGAAG': 'DDDDDDDDDDDDDDD<<:9' ,
#                                'AAAT': "DDDDDDDDDDAA@?====::876633220/-----*'" ,
#                                'AAATAAT': 'DDDDDDDDDDDDD@=;' ,
#                                'AAATT': 'DDDDDDDDDBAA>><:7' ,
#                                'AAC': 'DDDDDDDBA@<<<::988765553211/////,,,,,,,*)' ,
#                                'AACC': 'DDDDDDDDCB??=<<87633333332100/' ,
#                                'AACCCT': 'DDDDDDDDDDD>:755543210/.--++)&' ,
#                                'AAG': 'DDDDDDDBA@998444333221000//-,,++++++++++**)(' ,
#                                'AAGAG': 'DDDDDDDDDDCCA>;:887544320//.--,+' ,
#                                'AAGAGG': 'DDDDDDDDDDDCCA=;;8764321/' ,
#                                'AAGC': 'DDDDDDDDDCCA=:::987654' ,
#                                'AAGG': 'DDDDDDDDDCA@=<99988776665222000000////.---+++)' ,
#                                'AAGGAG': 'DDDDDDDDDDDAA@=;97544444332' ,
#                                'AAGGG': 'DDDDDDDDDCCCAA>=;:98654' ,
#                                'AAT': 'DDDDDDDDCB@?===<:97553322//..----,' ,
#                                'AATC': "DDDDDDDDDCCA>;::9999987666651//.-+*'" ,
#                                'AATG': 'DDDDDDDDDDDDDAAA???=<9777444443100/.....-,+)&' ,
#                                'AATGG': "DDDDDDDDDCA@<;888875554411100.........-+++++++****(((((('" ,
#                                'AATT': 'DDDDDDDDCCC@?>>==<;7542' ,
#                                'AC': "DDDDDDDB@>;87543210//...-,,,+++++****)))))))))))'" ,
#                                'ACACACACAG': 'DDDDDDDDDDDDDDDDDDD764322110/' ,
#                                'ACACACACAT': 'DDDDDDDDDDDDDDDDDDD55554' ,
#                                'ACACACACC': 'DDDDDDDDDDDDDDDDD5' ,
#                                'ACACACACGC': 'DDDDDDDDDDDDDDDDDDD5' ,
#                                'ACACACAG': 'DDDDDDDDDDDDDDD9877543' ,
#                                'ACACACAT': 'DDDDDDDDDDDDDDD98777643' ,
#                                'ACACACATATATAT': 'DDDDDDDDDDDDDDDDDDDDDDDDDDD420/.' ,
#                                'ACACACATGC': 'DDDDDDDDDDDDDDDDDDD::76420' ,
#                                'ACACACGC': 'DDDDDDDDDDDDDDD7742' ,
#                                'ACACACTC': 'DDDDDDDDDDDDDDD76' ,
#                                'ACACAG': 'DDDDDDDDDDDA@=<7654332' ,
#                                'ACACAT': 'DDDDDDDDDDD=;;:7666653222210/' ,
#                                'ACACATAT': 'DDDDDDDDDDDDDDD:9' ,
#                                'ACACATATAT': 'DDDDDDDDDDDDDDDDDDD9876' ,
#                                'ACACATATATAT': 'DDDDDDDDDDDDDDDDDDDDDDD444432' ,
#                                'ACACATGC': 'DDDDDDDDDDDDDDD<<;9' ,
#                                'ACACGC': 'DDDDDDDDDDD97731' ,
#                                'ACAG': "DDDDDDDDDC==;;87743331111110000//..-,,,+++*'&" ,
#                                'ACAGAG': 'DDDDDDDDDDDB===:998865432' ,
#                                'ACAGAGAG': 'DDDDDDDDDDDDDDD99976654/' ,
#                                'ACAGAGAGAG': 'DDDDDDDDDDDDDDDDDDD766543' ,
#                                'ACAT': 'DDDDDDDDBB>><::::9997775444421000/........-*)' ,
#                                'ACATAT': 'DDDDDDDDDDDA<<;;;:887775543200/.--,,,,+++++*' ,
#                                'ACATATAT': 'DDDDDDDDDDDDDDD=9' ,
#                                'ACATATATAT': 'DDDDDDDDDDDDDDDDDDD76' ,
#                                'ACC': 'DDDDDDDDDD?>>;:::776666211/' ,
#                                'ACCACCATC': 'DDDDDDDDDDDDDDDDD821+' ,
#                                'ACCATC': 'DDDDDDDDDDD@?==:8887554' ,
#                                'ACCCCC': 'DDDDDDDDDDDAAA=:' ,
#                                'ACCCCCC': 'DDDDDDDDDDDDD<9862' ,
#                                'ACCT': 'DDDDDDDDDDC?=<;::9865432110/' ,
#                                'ACGC': 'DDDDDDDA::6444443' ,
#                                'ACT': 'DDDDDDDDDAA@==:::8754220000/.---,,++****)' ,
#                                'ACTC': 'DDDDDDDDDDBB@=<<::876554' ,
#                                'ACTCC': 'DDDDDDDDDAA@@;855332' ,
#                                'ACTGCGG': 'DDDDDDDDDDDDD;' ,
#                                'AG': "DDDDDDCC??:96644211000.....-,++++*****)))))))'''''&%" ,
#                                'AGAGAGGG': 'DDDDDDDDDDDDDDD<<:8' ,
#                                'AGAGGG': 'DDDDDDDDDDDCA?<;::86555552' ,
#                                'AGAT': 'DDDDDDDDDB@?<;;866664320000///-,,,+++++++++**********)(' ,
#                                'AGATAGATAGATGAT': 'DDDDDDDDDDDDDDDDDDDDDDDDDDDDD00/' ,
#                                'AGATAGATGAT': 'DDDDDDDDDDDDDDDDDDDDD4310////..-,' ,
#                                'AGATAT': 'DDDDDDDDDDDB><<<<987666554432110//..--,' ,
#                                'AGC': "DDDDDDDDDDCA@@>=;;:7665531//------,,,**))'&$" ,
#                                'AGG': 'DDDDDDDDDCCA??><;9988865' ,
#                                'AGGC': 'DDDDDDDDDDDD@=<;84' ,
#                                'AGGG': 'DDDDDDDDDBB@?>;;;87' ,
#                                'AT': 'DDDDDDDB@=984310//.-,,,,,++++)(' ,
#                                'ATC': 'DDDDDDDDDDBB>>>=<9777765222111//--------,,,)' ,
#                                'ATCC': 'DDDDDDDDDDDAA====<;9777554222221//...---,,,,,,,,+*' ,
#                                'C': 'DDDDDB<62/..--,' ,
#                                'CCG': 'DDDDDDDDD@@@>><;98653210' ,
#                                'CG': 'DDDDDD?<7210//....-' ,
#                                'nonrepetitive': 'D' ,
#                                }

cdef dict indel_error_model = default_indel_error_model

###################################################################################################

cdef void my_free(void* thePointer):
    """
    Cython wrapper. Used for profiling.
    """
    free(thePointer)

###################################################################################################

cdef void* my_malloc(size_t theSize):
    """
    Cython wrapper. Used for profiling.
    """
    return malloc(theSize)

###################################################################################################

cdef void* my_calloc(size_t theSize1, size_t theSize2):
    """
    Cython wrapper. Used for profiling.
    """
    return calloc(theSize1, theSize2)

###################################################################################################

cdef void* my_realloc(void* thePointer, size_t newSize):
    """
    Cython wrapper. Used for profiling.
    """
    return realloc(thePointer, newSize)

###################################################################################################

cdef int computeOverlapOfReadAndHaplotype(int hapStart, int hapEnd, cAlignedRead* theRead):
    """
    Compute and return the number of bases by which a read overlaps the haplotype of interest.
    """
    cdef int readStart = theRead[0].pos
    cdef int readEnd = theRead[0].end
    cdef int overlapStart = max(hapStart, readStart)
    cdef int overlapEnd = min(hapEnd, readEnd)

    if overlapEnd > overlapStart:
        return overlapEnd - overlapStart
    else:
        return -1

###################################################################################################

cdef class Haplotype:
    """
    Class to encapsulate a single haplotype. This will store all the
    variants pertaining to this haplotype, as well as the reference
    sequence, all supporting reads, and the start and end positions of
    this haplotype in the reference.
    """
    def __init__(self, bytes refName, int startPos, int endPos, tuple variants, FastaFile refFile, int maxReadLength, int useIndelErrorModel, options):
        """
        Constructor. Takes a tuple of variants and a
        fasta file of the referene sequence.
        Variants are sorted on instantiation
        """
        self.refName = refName
        self.refFile = refFile
        self.variants = variants
        self.hash = -1
        self.indelErrorModel = indel_error_model
        self.useIndelErrorModel = useIndelErrorModel
        self.cLocalGapOpenQ = NULL
        self.haplotypeSequence = None
        self.startPos = max(0, startPos)
        self.endPos = min(endPos, self.refFile.refs[self.refName].SeqLength-1)
        self.maxReadLength = maxReadLength
        self.endBufferSize = min(2*maxReadLength, 200) # Cap the buffer size at a reasonable length
        self.verbosity = options.verbosity
        self.options = options
        self.lastIndividualIndex = -1

        cdef Variant v

        if len(variants) > 0:
            self.minVarPos = min([v.minRefPos for v in variants])
            self.maxVarPos = max([v.maxRefPos for v in variants])

            if self.minVarPos == self.maxVarPos:
                self.maxVarPos += 1
        else:
            self.minVarPos = self.startPos
            self.maxVarPos = self.endPos

        self.referenceSequence = self.refFile.getSequence(self.refName, self.startPos - self.endBufferSize, self.endPos + self.endBufferSize)

        if len(self.variants) == 0:
            self.haplotypeSequence = self.referenceSequence
        else:
            leftBuffer = self.refFile.getSequence(self.refName, self.startPos - self.endBufferSize, self.startPos)
            rightBuffer = self.refFile.getSequence(self.refName, self.endPos, self.endPos + self.endBufferSize)
            self.haplotypeSequence = leftBuffer + self.getMutatedSequence() + rightBuffer

        self.cHaplotypeSequence = self.haplotypeSequence
        self.hapLen = len(self.cHaplotypeSequence)

        #if self.referenceSequence == self.haplotypeSequence and len(self.variants) != 0:
        #    logger.error("Haplotype is broken. Var seq and ref seq are the same, and variants are %s" %(list(self.variants)))

        if self.hapLen > hash_size:
            logger.error("Haplotype with vars %s has len %s. Start is %s. End is %s. maxReadLen = %s" %(self.variants, self.hapLen, self.startPos, self.endPos, maxReadLength))
            logger.debug(self.haplotypeSequence)

        self.cHomopolQ = homopolq
        self.hapSequenceHash = NULL
        self.hapSequenceNextArray = NULL
        self.likelihoodCache = NULL
        self.lenCache = 0
        self.localGapOpen = NULL

    def __dealloc__(self):
        """
        Clean up cache.
        """
        if self.likelihoodCache != NULL:
            my_free(self.likelihoodCache)

        if self.hapSequenceHash != NULL:
            my_free(self.hapSequenceHash)

        if self.hapSequenceNextArray != NULL:
            my_free(self.hapSequenceNextArray)

        if self.localGapOpen != NULL:
            my_free(self.localGapOpen)

    def __copy__(self):
        """
        Make sure this never gets called for haplotypes.
        """
        raise StandardError, "Oh no! The bridge is gone!"

    def __richcmp__(Haplotype self, Haplotype other, int opCode):
        """
        Comparison function:

        Are two haplotypes equal? Only return true if the mutated
        sequences are exactly equal.
        """
        # <
        if opCode == 0:
            if self.refName < other.refName:
                return True
            elif self.refName == other.refName and self.startPos < other.startPos:
                return True
            elif self.refName == other.refName and self.startPos == other.startPos and self.haplotypeSequence < other.haplotypeSequence:
                return True
            else:
                return False
        # <=
        elif opCode == 1:
            if self.refName > other.refName:
                return False
            elif self.refName == other.refName and self.startPos > other.startPos:
                return False
            elif self.refName == other.refName and self.startPos == other.startPos and self.haplotypeSequence > other.haplotypeSequence:
                return False
            else:
                return True
        # >
        elif opCode == 4:
            if self.refName > other.refName:
                return True
            elif self.refName == other.refName and self.startPos > other.startPos:
                return True
            elif self.refName == other.refName and self.startPos == other.startPos and self.haplotypeSequence > other.haplotypeSequence:
                return True
            else:
                return False
        # >=
        elif opCode == 5:
            if self.refName < other.refName:
                return False
            elif self.refName == other.refName and self.startPos < other.startPos:
                return False
            elif self.refName == other.refName and self.startPos == other.startPos and self.haplotypeSequence < other.haplotypeSequence:
                return False
            else:
                return True
        # ==
        if opCode == 2:
            if self.refName != other.refName:
                return False
            elif self.startPos != other.startPos:
                return False
            else:
                thisSeq = self.haplotypeSequence
                otherSeq = other.haplotypeSequence
                return thisSeq == otherSeq

        # !=
        elif opCode == 3:
            if self.refName != other.refName:
                return True
            elif self.startPos != other.startPos:
                return True
            else:
                thisSeq = self.haplotypeSequence
                otherSeq = other.haplotypeSequence
                return thisSeq != otherSeq
        else:
            raise StandardError, "Op code %s not implemented in haplotype__richcmp__()" %(opCode)

    def __hash__(self):
        """
        Implementing this function allows haplotypes to be hashed, and so stored in
        a set or dictionary. The supporting reads are not included in the hashing, as
        we want two haplotypes to give the same hash id if they have the same positions
        and sequences.
        """
        if self.hash == -1:
            self.hash = hash((self.refName, self.startPos, self.endPos, self.haplotypeSequence))

        return self.hash

    cdef double* alignReads(self, int individualIndex, cAlignedRead** start, cAlignedRead** end, cAlignedRead** badReadsStart, cAlignedRead** badReadsEnd, cAlignedRead** brokenReadsStart, cAlignedRead** brokenReadsEnd, int useMapQualCap, int printAlignments):
        """
        """
        cdef int readIndex = 0
        cdef double score = 0.0
        cdef int nReads = end - start
        cdef int nBadReads = badReadsEnd - badReadsStart
        cdef int nBrokenReads = brokenReadsEnd - brokenReadsStart
        cdef int totalReads = nReads + nBadReads + nBrokenReads
        cdef int readOverlap = 0
        cdef int hapLen = self.endPos - self.startPos
        cdef int readLen = 0
        cdef double* temp = NULL

        # Either first time, or new individual
        if individualIndex != self.lastIndividualIndex:

            if self.likelihoodCache == NULL:
                self.likelihoodCache = <double*>(my_malloc((totalReads+1)*sizeof(double)))
                self.lenCache = totalReads

                if self.likelihoodCache == NULL:
                    logger.error("Could not allocate haplotype cache")
                    raise StandardError, "Out of memory in cHaplotype.alignReads"
            else:
                if totalReads >= self.lenCache:
                    temp = <double*>realloc(self.likelihoodCache, 2*totalReads*sizeof(double))

                    if temp == NULL:
                        logger.error("Could not reallocate haplotype cache")
                        raise StandardError, "Out of memory in cHaplotype.alignReads"

                    self.likelihoodCache = temp
                    self.lenCache = 2*totalReads

            self.lastIndividualIndex = individualIndex

            if printAlignments:
                logger.debug("")
                logger.debug("#########################################################################")
                logger.debug("Logging alignments for haplotype %s and sample %s" %(self, individualIndex))
                logger.debug("Haplotype sequence is followed by alignments of all reads (good then bad then broken mates)")
                logger.debug("#########################################################################")
                logger.debug(self.getMutatedSequence())
                logger.debug("")

            while start != end:

                readOverlap = computeOverlapOfReadAndHaplotype(self.startPos, self.endPos, start[0])

                if Read_IsQCFail(start[0]) or readOverlap < hash_nucs:
                    self.likelihoodCache[readIndex] = 0
                else:
                    score = alignReadToHaplotype(start[0], self, useMapQualCap, printAlignments)
                    self.likelihoodCache[readIndex] = score

                start += 1
                readIndex += 1

            while badReadsStart != badReadsEnd:

                readOverlap = computeOverlapOfReadAndHaplotype(self.startPos, self.endPos, badReadsStart[0])

                if Read_IsQCFail(badReadsStart[0]) or readOverlap < hash_nucs:
                    self.likelihoodCache[readIndex] = 0
                else:
                    score = alignReadToHaplotype(badReadsStart[0], self, useMapQualCap, printAlignments)
                    self.likelihoodCache[readIndex] = score

                badReadsStart += 1
                readIndex += 1

            # It doesn't make sense to check overlap for the broken mates, as their mapping positions don't make
            # sense in this context.
            while brokenReadsStart != brokenReadsEnd:
                score = alignReadToHaplotype(brokenReadsStart[0], self, useMapQualCap, printAlignments)
                self.likelihoodCache[readIndex] = score
                brokenReadsStart += 1
                readIndex += 1

            self.likelihoodCache[readIndex] = 999 # End marker

        return self.likelihoodCache

    cdef inline double alignSingleRead(self, cAlignedRead* theRead, int useMapQualCap, int printAlignments):
        """
        Returns the alignment score for a single read. If 'useMapQualCap' is True, then read likelihood
        is capped using the mapping quality of the read. Otherwise it is capped at 1e-300.
        """
        return alignReadToHaplotype(theRead, self, useMapQualCap, printAlignments)

    cdef char* getReferenceSequence(self, prefix = 0):
        """
        Return the refernece sequence for the region covered by this haplotype. pretty shows where the variants are.
        """
        if prefix == 0 and self.referenceSequence != None:
            return self.referenceSequence

        seqMax = self.refFile.refs[self.refName].SeqLength - 1
        self.referenceSequence = self.refFile.getSequence(self.refName, max(0, self.startPos - prefix), min(self.endPos + prefix, seqMax))
        return self.referenceSequence

    cdef char* getMutatedSequence(self):
        """
        Return the reference sequence mutated with all the variants being
        considered in this haplotype.

        Just to remind ourselves: SNPs are reported at the index of the changed
        base (0-indexed internally, and 1-indexed in VCF). Insertions are reported
        such that the index is that of the last reference base before the insertion.
        Deletions are reported such that the index is that of the first deleted base.
        """
        cdef Variant v
        cdef Variant firstVar

        if self.haplotypeSequence is None:

            currentPos = self.startPos

            # Get sequence up to one base before the first variant
            firstVar = self.variants[0]
            bitsOfMutatedSeq = [self.refFile.getSequence(self.refName, currentPos, firstVar.refPos)]
            currentPos = firstVar.refPos

            for v in self.variants:

                # Move up to one base before the next variant, if we're not already there.
                if v.refPos > currentPos:
                    bitsOfMutatedSeq.append(self.refFile.getSequence(self.refName, currentPos, v.refPos))
                    currentPos = v.refPos

                # SNP/Mult-SNP/Complex
                if v.nAdded == v.nRemoved:
                    bitsOfMutatedSeq.append(v.added)
                    currentPos += v.nRemoved

                # Arbitrary length-changing sequence replacement
                else:
                    if v.refPos == currentPos:
                        bitsOfMutatedSeq.append(self.refFile.getCharacter(self.refName, v.refPos))
                        currentPos += 1

                    currentPos += v.nRemoved
                    bitsOfMutatedSeq.append(v.added)

            # Is this ok when currentPos == endPos?
            if currentPos > self.endPos:
                logger.error("cpos = %s end pos = %s. Variants are %s" %(currentPos, self.endPos, self.variants))

            if currentPos < self.endPos:
                bitsOfMutatedSeq.append(self.refFile.getSequence(self.refName, currentPos, self.endPos))

            self.haplotypeSequence = bytes(''.join(bitsOfMutatedSeq))

        return self.haplotypeSequence

    cdef list homopolymerLengths( self ):
        """
        Return a list of homopolymer lengths for the sequence surrounding
        each variant in this haplotype.
        """
        cdef tuple vsf = self.variants
        if len( vsf ) == 0:
            return []
        else:
            return [ self.homopolymerLengthForOneVariant( v ) for v in vsf ]

    cdef int homopolymerLengthForOneVariant(self, Variant variant):
        """
        Calculate and return the length of the largest homopolymer
        touching this variant. Compute homopolymer lengths on the
        left and right, and return the largest.
        """
        varChrom = variant.refName
        varPos = variant.refPos

        leftRefSeq = self.refFile.getSequence(varChrom, varPos-20, varPos)
        rightRefSeq = self.refFile.getSequence(varChrom, varPos+1, varPos + 21)

        if len(leftRefSeq) == 0 or len(rightRefSeq) == 0:
            return 0

        leftHpSize = 0
        rightHpSize = 0

        firstLeftChar = leftRefSeq[-1]
        firstRightChar = rightRefSeq[0]

        for char in reversed(leftRefSeq):
            if char == firstLeftChar:
                leftHpSize += 1
            else:
                break

        for char in rightRefSeq:
            if char == firstRightChar:
                rightHpSize += 1
            else:
                break

        if firstLeftChar != firstRightChar:
            return max(leftHpSize, rightHpSize)
        else:
            return leftHpSize + rightHpSize

    cdef bytes getSequenceContext(self, Variant variant):
        """
        Return the sequence surrounding this variant's position.
        """
        varChrom = variant.refName
        varPos = variant.refPos
        return self.refFile.getSequence(varChrom, varPos-10, varPos + 11)

    cdef dict vcfINFO(self):
        """
        Information to go in the vcf file INFO field - a two level dictionary of variant - field - value
        This can be augmented at the individual/population level with more information, or some of it can be
        removed before printing the output

        Include
        HP = Honopolymer tract length
        NR = Number of supporting reads
        CD = Coverage depth
        """
        cdef Variant variant
        cdef dict INFO = {}

        homopolymerLengths = self.homopolymerLengths()
        vsf = self.variants

        for varIndex, variant in enumerate( self.variants ):
            HP = homopolymerLengths[varIndex]
            SC = self.getSequenceContext(variant)

            INFO[vsf[varIndex]] = {'HP':[HP], 'SC':[SC]}

        return INFO

    def __str__(self):
        """
        Generate a string representation of the haplotype. This is useful for
        debugging.
        """
        if len(self.variants) == 0:
            return '  Haplotype(*Reference*) %s:%s-%s' %(self.refName, self.startPos, self.endPos)

        vars = [str(v) for v in self.variants]
        string = "  Haplotype(" + ",".join(vars) + ") %s:%s-%s" %(self.refName, self.startPos, self.endPos)

        return string

    def __repr__(self):
        """
        The representation function. Called when printing the screen.
        """
        return self.__str__()

    cdef void annotateWithGapOpen(self):
        """
        Annotate this haplotype with a context-specific gap open penalty, using either the
        homopolymer model or the more complex indel model, as specified by the self.useIndelErrorModel
        flag.
        """
        # Only do this one per haplotype
        if self.localGapOpen != NULL:
            return

        cdef int index = self.hapLen
        cdef char* seq = self.cHaplotypeSequence
        cdef char* errorModel = NULL
        cdef int homopol = -1
        cdef int homopollen = 0

        self.localGapOpen = <short*>malloc(self.hapLen*sizeof(short))

        if self.useIndelErrorModel and self.cLocalGapOpenQ == NULL:
            hapSeqForAnnotating = "".join([c if c != "N" else random.choice(["A", "C", "T", "G"]) for c in self.haplotypeSequence])
            self.localGapOpenQ = cerrormodel.annotate_sequence(hapSeqForAnnotating, self.indelErrorModel, 0)
            self.cLocalGapOpenQ = self.localGapOpenQ

        if self.useIndelErrorModel:
            errorModel = self.cLocalGapOpenQ
        else:
            errorModel = self.cHomopolQ

        if not self.useIndelErrorModel:

            homopol = -1
            homopollen = 0

            while index > 0:

                index -= 1

                if seq[index] == homopol:
                    homopollen += <int>(not(not(errorModel[homopollen+1])))
                else:
                    homopollen = 0

                self.localGapOpen[index] = 4*(<int>(errorModel[homopollen]) - (<int>'!'))
                homopol = seq[index];

                if homopol == 'N':
                    homopol = 0

        else:
            while index > 0:
                index -= 1
                self.localGapOpen[index] = (<int>errorModel[index])*4

###################################################################################################

cdef double alignReadToHaplotype(cAlignedRead* read, Haplotype hap, int useMapQualCap, int printAlignments):
    """
    This is the basic, banded-alignment routine that forms the heart of Platypus. This function decides where to anchor the
    read sequence to the specified haplotype, and calls the fastAlignmentRoutine function, which performs a banded alignment.

    If we don't anchor the read sequence to the correct part of the haplotype, then all the results, particularly for indels,
    will be rubbish.
    """
    cdef char* hapSeq = hap.cHaplotypeSequence
    cdef char* aln1 = NULL
    cdef char* aln2 = NULL
    cdef char* errorModel = NULL
    cdef int hapStart = hap.startPos - hap.endBufferSize
    cdef int gapExtend = 3
    cdef int nucprior = 2
    cdef int hapLen = hap.hapLen
    cdef int useHomopolQ = True

    if hap.useIndelErrorModel and hap.cLocalGapOpenQ == NULL:
        hapSeqForAnnotating = "".join([c if c != "N" else random.choice(["A", "C", "T", "G"]) for c in hap.haplotypeSequence])
        hap.localGapOpenQ = cerrormodel.annotate_sequence(hapSeqForAnnotating, hap.indelErrorModel, 0)
        hap.cLocalGapOpenQ = hap.localGapOpenQ

    if hap.useIndelErrorModel:
        errorModel = hap.cLocalGapOpenQ
        useHomopolQ = False
    else:
        errorModel = hap.cHomopolQ
        useHomopolQ = True

    cdef char* readSeq = read.seq
    cdef char* readQuals = read.qual

    cdef int readStart = read.pos
    cdef int readLen = read.rlen
    cdef int mapQual = read.mapq

    cdef int lenOfHapSeqToTest = readLen + 15
    cdef int hapPos = readStart - hapStart
    cdef int alignScore = 0

    cdef double probMapWrong = mLTOT*read.mapq  # A log value
    cdef double probMapRight = log(1.0 - exp(mLTOT*read.mapq)) # A log value

    # Arbitrary cap
    cdef double likelihoodCap = 0.0

    if useMapQualCap == True:
        likelihoodCap = probMapWrong
    else:
        likelihoodCap = -300 # Arbitrary cap close to min float value

    # Make sure read hash exists when needed
    if read.hash == NULL:
        hashReadForMapping(read)

    # Make sure haplotype hash exists when needed
    if hap.hapSequenceHash == NULL:
        hash_sequence_multihit(hap.cHaplotypeSequence, hap.hapLen, &hap.hapSequenceHash, &hap.hapSequenceNextArray)

    cdef int mapPos = -1

    if printAlignments:
        aln1 = <char*>(calloc(2*(readLen + hapLen), sizeof(char)))
        aln2 = <char*>(calloc(2*(readLen + hapLen), sizeof(char)))

        for i in range(2*(readLen + hapLen)):
            aln1[i] = 'N'
            aln2[i] = 'N'

    alignScore = mapAndAlignReadToHaplotype(readSeq, readQuals, readStart, hapStart, readLen, hapLen, hap.hapSequenceHash, hap.hapSequenceNextArray, read.hash, hapSeq, gapExtend, nucprior, errorModel, useHomopolQ, aln1, aln2, &mapPos)

    if printAlignments:
        logAlignmentOfReadToHaplotype(hap, readSeq, readQuals, aln1, aln2, mapPos, alignScore, hap.options)
        free(aln1)
        free(aln2)

    return max(mLTOT*alignScore + probMapRight, likelihoodCap)

###################################################################################################

cdef void logAlignmentOfReadToHaplotype(Haplotype hap, char* readSeq, char* readQuals, char* aln1, char* aln2, int mapPos, int alignScore, options):
    """
    Does exactly what it says on the tin: uses the logger to output a string representation of
    the final alignment of this read to this haplotype.
    """
    cdef list alignmentChars = [" "] * mapPos
    cdef str hapSeq = str(hap.haplotypeSequence)
    cdef str aln1Str = str(aln1)
    cdef str aln2Str = str(aln2)
    cdef int minQual = options.minBaseQual
    cdef int qualIndex = 0
    cdef int hapIndex = mapPos
    cdef int hadInsertion = False

    for theChar1,theChar2 in zip(aln1Str,aln2Str):

        if hapIndex >= len(hapSeq):
            break

        # Gap in read: deletion of ref bases.
        if theChar2 == "-":
            alignmentChars.append("-")
            hapIndex += 1

        # Gap in hap: deletion of ref bases.
        elif theChar1 == "-":
            #alignmentChars.append("-")
            #hapOffset -= 1
            qualIndex += 1
            hadInsertion = True

        elif theChar2 == hapSeq[hapIndex]:

            if hadInsertion:
                alignmentChars.append("!")
            else:
                alignmentChars.append(".")

            hadInsertion = False

            qualIndex += 1
            hapIndex += 1

        else:
            newChar = None

            if readQuals[qualIndex] > minQual:
                newChar = theChar2.lower()
            else:
                newChar = "n"

            if hadInsertion:
                newChar = newChar.upper()

            alignmentChars.append( newChar )
            hadInsertion = False

            qualIndex += 1
            hapIndex += 1

    logger.debug("".join(alignmentChars) + " Score = %s" %(alignScore))

###################################################################################################

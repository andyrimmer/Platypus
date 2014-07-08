seq = "TATTTGCATGCGCTTTCGAGCTGTTGAAGAGACGTGTATTGGAATAAGTAATCACATAAGTGTTAGTAACTTATTTAAATACGTATAGAGTCGCCTATTTGCCTAGCCTTTTGGTTCTCAGATTTTTTAATTATTACATTGCTATAAGGGTGTAACTGTGTGATAGCCAAAATTTTAAGCTGCAAATGGTTTGTAAATATGATATATTACAAGCTTCATGAAAATCGGTTTATGACTGATCCGCGATTACGTTGAAAGGCGACTGGCAGAGATACTTTTGTTCAGATGTTTTTTCAGGTAGCGATTCCAATGAATAGGTAAAATACCTTGCAAGTTTTGTTGTTGTCGTTGGAGGAAATGTGGATGTGGTTGTTATTGTTGA"

prepost = "AAAAAAAAAAAAAAAAAA"
read1 = seq[100:150]
read2 = seq[0:50]
read3 = prepost+seq[0:50]
read4 = seq[-50:]
read5 = seq[-40:] + prepost

import calign

h = calign.hash_sequence(seq)
print "Testing..."
assert calign.map_read(read1, len(seq), h) == 100
assert calign.map_read(read2, len(seq), h) == 0
assert calign.map_read(read3, len(seq), h) == -len(prepost)
assert calign.map_read(read4, len(seq), h) == len(seq)-50
assert calign.map_read(read5, len(seq), h) == len(seq)-40
assert calign.map_cortex(read2 + read1, seq, len(seq), h, 10) == (0,0,50,50,100,50)
assert calign.map_cortex(read2 + prepost + read1, seq, len(seq), h, 10) == (0,0,51,50+len(prepost),100,50)
# 51 because first A of prepost matches
assert calign.map_cortex(prepost + read2 + prepost, seq, len(seq), h, 10) == (len(prepost),0,51)
# 51 because first A of prepost matches
assert calign.map_cortex(prepost, seq, len(seq), h, 10) == None
print "Test Complete"

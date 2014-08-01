###################################################################################################

cdef extern from "stdlib.h":
    void free(void *)
    void *malloc(size_t)
    void *realloc(void *,size_t)
    void *calloc(size_t,size_t)
    void *memset(void *buffer, int ch, size_t count )
    void qsort(void*, size_t, size_t, int(*)(void*, void*))

###################################################################################################

# define the minimum tandem length defined in the model.
cdef int MINIMUM_TANDEM_LENGTH = 4
cdef dict repunit_dict = None

cdef extern from "tandem.h":
    void annotate( char* sequence, char* sizes, char* displacements, int length)
    int approximate_indel_rate(int size, int displacement)

###################################################################################################

cdef tuple calculate_size_and_displacement( bytes sequence, char annotate_all ):
    cdef char* csequence = sequence
    cdef int length = len(sequence)
    cdef char* csizes = <char*>calloc(length+1, sizeof(char))
    cdef char* cdisplacements = <char*>calloc(length+1, sizeof(char))
    if annotate_all:
        annotate( csequence, csizes, cdisplacements, -length )
    else:
        annotate( csequence, csizes, cdisplacements, length )
    cdef bytes sizes = csizes
    cdef bytes displacements = cdisplacements
    free (csizes)
    free (cdisplacements)
    return (sizes, displacements)

###################################################################################################

cdef bytes normalize_repunit_fast( bytes unit ):
    global repunit_dict
    cdef int length = len(unit)
    cdef char* cunit = unit
    cdef int i
    cdef int both
    for i in range(length):
        cunit[i] &= 0xDF
        if cunit[i] == 'N':
            return <bytes>("N" * length)
    if length == 1:
        if cunit[0] <= 'C': return unit
        if cunit[0] == 'G': return bytes('c')
        return bytes('a')
    elif length == 2:
        both = <int>cunit[0]+<int>cunit[1]
        if both == 138: return bytes('CG')   # c+g
        elif both == 149: return bytes('AT') # a+t
        elif both == 132: return bytes('AC') # a+c 
        elif both == 155: return bytes('ac') # g+t
        elif both == 136: return bytes('AG') # a+g
        return bytes('ag') # t+c
    cdef bytes n1,n2,n3,n4,n5,n6,n7,n8
    if repunit_dict == None:
        repunit_dict = {}
        for n1 in bytes('ACGT'):
            for n2 in bytes('ACGT'):
                for n3 in bytes('ACGT'):
                    repunit_dict[n1+n2+n3] = normalize_repunit(n1+n2+n3)
                    for n4 in bytes('ACGT'):
                        repunit_dict[n1+n2+n3+n4] = normalize_repunit(n1+n2+n3+n4)
                        for n5 in bytes('ACGT'):
                            repunit_dict[n1+n2+n3+n4+n5] = normalize_repunit(n1+n2+n3+n4+n5)
                            for n6 in bytes('ACGT'):
                                repunit_dict[n1+n2+n3+n4+n5+n6] = normalize_repunit(n1+n2+n3+n4+n5+n6)
                                for n7 in bytes('ACGT'):
                                    repunit_dict[n1+n2+n3+n4+n5+n6+n7] = normalize_repunit(n1+n2+n3+n4+n5+n6+n7)
                                    for n8 in bytes('ACGT'):
                                        repunit_dict[n1+n2+n3+n4+n5+n6+n7+n8] = normalize_repunit(n1+n2+n3+n4+n5+n6+n7+n8)
    if length <= 8:
        return repunit_dict[cunit]
    return normalize_repunit(cunit)

###################################################################################################

cdef bytes normalize_repunit( bytes unit ):
    cdef int length = len(unit)
    cdef list unit2list = [{'A':'T','T':'A','C':'G','G':'C'}.get(c,'N') for c in unit]
    unit2list.reverse()
    cdef bytes unit2 = bytes(''.join(unit2list))
    unit2 += unit2
    unit += unit
    cdef bytes normunit = sorted([ unit[i:i+length] for i in range(length) ] + [ unit2[i:i+length]+'-' for i in range(length) ])[0]
    if normunit[-1] == '-':
        normunit = normunit[:-1].lower()
    return normunit

###################################################################################################

cpdef list get_repeats( bytes sequence, int min_length, long pos ):
    """ Returns list of non-overlapping (position, size, tandemunit); 
        tandemunit is lower case for units on the reverse strand """

    cdef char* seq = sequence
    cdef int seqlen = len(sequence)
    cdef bytes sizes
    cdef bytes displacements
    cdef char* csizes
    cdef char* cdisplacements
    cdef list repeats
    cdef int idx, size, displacement
    cdef bytes tandemunit

    (sizes, displacements) = calculate_size_and_displacement( sequence, True )
    csizes = sizes
    cdisplacements = displacements

    repeats = []

    for idx in range(seqlen):

        size = <int>csizes[idx]
        if size >= min_length:

            # only add if not overlapping previously entered repeat
            displacement = <int>cdisplacements[idx]
            if len(repeats)==0 or repeats[-1][1] != size or len(repeats[-1][2]) != displacement or repeats[-1][0] + size < pos+idx:
                tandemunit = normalize_repunit_fast(sequence[ idx : idx+displacement ])
                if tandemunit[0] != 'N':
                    repeats.append( (pos+idx, size, tandemunit) )

    return repeats

###################################################################################################

cpdef int approx_indel_error_rate( int size, int displacement ):
    """ Returns estimated indel error rate; larger is more """
    return approximate_indel_rate(size, displacement)

###################################################################################################

cdef tuple try_microsat( char* seq, int seqlen, int start, int replen, int direction ):
    cdef int minpos = start
    cdef int maxpos = start
    cdef int length
    cdef bytes character

    while minpos >= 0 and 0 <= minpos+replen*direction < seqlen and seq[minpos] == seq[minpos + replen*direction] and seq[minpos] != 'N':
        minpos -= 1

    while maxpos < seqlen and maxpos+replen*direction < seqlen and seq[maxpos] == seq[maxpos + replen*direction] and seq[maxpos] != 'N':
        maxpos += 1

    length = maxpos - minpos + replen - 1

    if length < 2*replen:
        character = seq[start]
        return 1, character, start

    return length, seq[minpos+1:minpos+replen+1], min(minpos+1,minpos+1+replen*direction)

###################################################################################################

cdef tuple microsattelite( char* seq, int seqlen, int start, int tandem ):
    cdef tuple maxlen = (-1,"",-1)
    cdef tuple result
    for m in range(1,tandem+1):
        result = try_microsat(seq, seqlen, start, m, -1)

        if result[0] > maxlen[0]: 
            maxlen = result
        result = try_microsat(seq, seqlen, start, m, 1)

        if result[0] > maxlen[0]: 
            maxlen = result
    return maxlen

###################################################################################################

def get_annotation(seq, pos, tandem):
    t, repunitT, repstart = microsattelite( seq, len(seq), pos, tandem )
    return t,normalize_repunit(repunitT),repstart

###################################################################################################

def add_tandem(int pos, int tandemlen, bytes tandemunit, list indelq, dict indel_q_data, int output_base=0):
    """
    Adds gap opening penalties in string indelq, for the tandem repeat described by output.
    indel_q_data contains a dictionary of error models.  These have the following two forms:
    {'AAC' : '#######',       # gap opening penalties, phred 33-based, for AAC repeats of length 1,2,..,7
     3 : '######'}            # gap opening penalties for any other triplet repeats of length 1,2,...,6
    Nonrepetitive rates are encoded in the first character of the entry for '1'.
    """
    tandemunit = tandemunit.upper()

    if pos == -1:
        return

    cdef bytes model = None
    cdef char* cModel = NULL
    cdef int qdata = 99

    if tandemunit in indel_q_data:
        model = indel_q_data[tandemunit]
        cModel = model
        qdata = cModel[ min(tandemlen-1, len(cModel)-1) ]

    # use generic model when specific model not available, or
    # when extrapolated specific rates are used, and generic model is available, and predicts higher indel rate
    if len(tandemunit) in indel_q_data and ((model is None) or (tandemlen > len(cModel))):
        model = indel_q_data[ len(tandemunit) ]
        cModel = model
        qdata = min( qdata, cModel[ min(tandemlen-1, len(cModel)-1) ] )

    if qdata == 99:
        return

    cdef int q = qdata - 33 + output_base

    for idx in range(pos, pos+tandemlen):
        indelq[idx] = chr( min(q, ord(indelq[idx]) ) )

###################################################################################################


def add_tandem_read(int pos, int tandemlen, bytes tandemunit, list indelq, dict indel_q_data, int output_base=0):
    """
    As add_tandem above, but calculates upper bounds for read base qualities due to polymerase slippage.
    The read is assumed to be oriented in the forward direction.
    Only quality scores overlaid by the tandem repeat are altered; these need to be extended rightward by the caller.
    """
    tandemunit = tandemunit.upper()

    if pos == -1:
        return

    cdef bytes model = None
    cdef bytes genericModel = None
    cdef char* cModel = NULL
    cdef char* cGenericModel = NULL
    cdef int qdata = 99
    cdef int modelidx
    cdef int q

    if tandemunit in indel_q_data:
        model = indel_q_data[tandemunit]
        cModel = model

    if len(tandemunit) in indel_q_data:
        genericModel = indel_q_data[len(tandemunit)]
        cGenericModel = genericModel

    for idx in range(pos, pos+tandemlen):

        qdata = 99
        modelidx = idx - pos

        if model is not None:
            qdata = cModel[ min(modelidx, len(cModel)-1) ]

        # use generic model when specific model not available, or
        # when extrapolated specific rates are used, and generic model is available, and predicts higher indel rate
        if genericModel is not None and ((model is None) or (modelidx > len(cModel))):
            qdata = min( qdata, cGenericModel[ min(modelidx, len(cGenericModel)-1) ] )

        if qdata == 99:
            return

        indelq[idx] = chr( min(qdata - 33 + output_base, ord(indelq[idx])) )

###################################################################################################


cdef bytes annotate_sequence_slow(sequence, dict indel_q_data, int output_base, int tandem):
    """
    Annotates a sequence with the local indel probability (gap opening penalty),
    using an error model described in indel_q_data
    """
    pos = 0
    indelq = [chr(ord(indel_q_data[1][0])-ord('!')+output_base)] * len(sequence)
    oldoutput = [-1,-1,-1]

    while pos < len(sequence):

        tandemlen, tandemunit, repstart = get_annotation(sequence, pos, tandem)
        output = [repstart, tandemlen, tandemunit]

        # concatenate or output
        if tandemlen >= 2:
            if oldoutput[0] + oldoutput[1] >= output[0] and oldoutput[2] == output[2]:
                oldoutput[1] = output[0] + output[1] - oldoutput[0]
            else:
                if oldoutput[0] != -1:
                    add_tandem(oldoutput[0], oldoutput[1], oldoutput[2], indelq, indel_q_data, output_base = output_base)
                oldoutput = output
            pos = repstart + tandemlen
        else:
            pos += 1

    if oldoutput[0] != -1:
        add_tandem(oldoutput[0], oldoutput[1], oldoutput[2], indelq, indel_q_data, output_base = output_base)

    return bytes(''.join(indelq))

###################################################################################################

cdef bytes annotate_sequence(sequence, dict indel_q_data, int output_base):
    """
    Annotates a sequence with the local indel probability (gap opening penalty),
    using an error model described in indel_q_data.
    """
    cdef list indelq = [chr(ord(indel_q_data[1][0])-ord('!')+output_base)] * len(sequence)
    cdef bytes sizes = None
    cdef bytes displacements = None

    (sizes, displacements) = calculate_size_and_displacement( sequence, False )

    cdef char* csizes = sizes
    cdef char* cdisplacements = displacements
    cdef int length = len(sequence)
    cdef int pos = 0
    cdef int tandemunitlen = -1
    cdef int tandemlen = -1
    cdef int oldoutputpos = -1
    cdef int oldoutputlen = -1
    cdef bytes tandemunit = None
    cdef bytes oldoutputunit = None

    while pos < length:

        tandemunitlen = cdisplacements[pos]
        tandemlen = csizes[pos]
        tandemunit = normalize_repunit_fast(sequence[pos:pos+tandemunitlen])

        # concatenate or output
        if tandemlen >= 2 and tandemunit.count('N') == 0:
            if oldoutputpos + oldoutputlen >= pos and oldoutputunit == tandemunit:
                oldoutputlen = pos + tandemlen - oldoutputpos
            else:
                if oldoutputpos != -1 and oldoutputlen >= MINIMUM_TANDEM_LENGTH:
                    add_tandem( oldoutputpos, oldoutputlen, oldoutputunit, indelq, indel_q_data, output_base = output_base)

                oldoutputpos = pos
                oldoutputlen = tandemlen
                oldoutputunit = tandemunit

        # next position
        pos += 1

    # add last tandem
    if oldoutputpos != -1:
        add_tandem( oldoutputpos, oldoutputlen, oldoutputunit, indelq, indel_q_data, output_base = output_base)

    # return result
    return bytes(''.join(indelq))


###################################################################################################

cdef bytes annotate_sequence_read(sequence, dict indel_q_data, int output_base, int reverse):
    """
    Annotates a sequence with the local indel probability (gap opening penalty),
    using an error model described in indel_q_data.
    """
    if reverse:
        sequence = ''.join(reversed(sequence))

    cdef list indelq = [chr(ord(indel_q_data[1][0])-ord('!')+output_base)] * len(sequence)
    cdef bytes sizes = None
    cdef bytes displacements = None

    (sizes, displacements) = calculate_size_and_displacement( sequence, False )

    cdef char* csizes = sizes
    cdef char* cdisplacements = displacements
    cdef int length = len(sequence)
    cdef int pos = 0
    cdef int tandemunitlen = -1
    cdef int tandemlen = -1
    cdef int oldoutputpos = -1
    cdef int oldoutputlen = -1
    cdef bytes tandemunit = None
    cdef bytes oldoutputunit = None
    cdef char q

    while pos < length:

        tandemunitlen = cdisplacements[pos]
        tandemlen = csizes[pos]
        tandemunit = normalize_repunit_fast(sequence[pos:pos+tandemunitlen])

        # concatenate or output
        if tandemlen >= 2 and tandemunit.count('N') == 0:
            if oldoutputpos + oldoutputlen >= pos and oldoutputunit == tandemunit:
                oldoutputlen = pos + tandemlen - oldoutputpos
            else:
                if oldoutputpos != -1 and oldoutputlen >= MINIMUM_TANDEM_LENGTH:
                    add_tandem_read( oldoutputpos, oldoutputlen, oldoutputunit, indelq, indel_q_data, output_base = output_base)

                oldoutputpos = pos
                oldoutputlen = tandemlen
                oldoutputunit = tandemunit

        # next position
        pos += 1

    # add last tandem
    if oldoutputpos != -1:
        add_tandem_read( oldoutputpos, oldoutputlen, oldoutputunit, indelq, indel_q_data, output_base = output_base)

    # post-process
    q = ord(indelq[0])
    for idx in range(len(indelq)):
        indelq[idx] = chr( min(q, ord(indelq[idx])))
        q = ord(indelq[idx])

    # return result
    if reverse:
        return bytes(''.join(reversed(indelq)))
    else:
        return bytes(''.join(indelq))


###################################################################################################

cpdef testAnnotate():
    seq1 = "TATTTGAAAAAAAAAAACATGCGCTTTCGAGCTGTTGAAGAGACGTGTATTGGAATAAGTAATCACATAAGTGTTAGTAACTTATTTAAATACGTATAGAGTCGCCTATTTGCCTAGCCTTTTGGTTCTCAGATTTTTTAATTATTACATTGCTATAAGGGTGTAACTGTGTGATAGCCAAAATTTTAAGCTGCAAATGGTTTGTAAATATGATATATTACAAGCTTCATGAAAATCGGTTTATGACTGATCCGCGATTACGTTGAAAGGCGACTGGCAGAGATACTTTTGTTCAGATGTTTTTTCAGGTAGCGATTCCAATGAATAGGTAAAATACCTTGCAAGTTTTGTTGTTGTCGTTGGAGGAAATGTGGATGTGGTTGTTATTGTTGA"

    import time

    t = time.clock()
    model = {1:'SSI?5+#',
             'AG':'SS#'}
    indelq = annotate_sequence_slow(seq1, model, ord('!'), 24)
    indelq_fast = annotate_sequence(seq1, model, ord('!'))
    asr_fw = annotate_sequence_read(seq1, model, ord('!'), False)
    asr_bw = annotate_sequence_read(seq1, model, ord('!'), True)
    print seq1
    print indelq
    print indelq_fast
    print asr_fw
    print asr_bw

    # test speed
    seq2 = seq1*10
    for i in range(10000):
        indelq = annotate_sequence(seq1, model, ord('!'))

    print time.clock()-t



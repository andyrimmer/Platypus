
cdef bytes annotate_sequence(sequence, dict indel_q_data, int output_base)
cdef bytes annotate_sequence_read(sequence, dict indel_q_data, int output_base, int reverse)
cdef bytes annotate_sequence_slow(sequence, dict indel_q_data, int output_base, int tandem)
cdef tuple calculate_size_and_displacement( bytes sequence, char annotate_all )
cpdef list get_repeats( bytes sequence, int min_length, long pos )
cpdef int approx_indel_error_rate( int size, int displacement )

cpdef testAnnotate()

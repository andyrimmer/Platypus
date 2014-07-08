"""
Here I'll dump all the utility code reqiured for creating useful summary
statistics. This will be used to inform hard-ref calls and similar
"""

###################################################################################################

class LocusSummary(object):
    """
    Stores summary statistics of sequence coverage at a
    specific locus.
    """
    def __init__(self, chrom, pos, refBase):
        """
        Constructor.
        """
        self.chrom = chrom
        self.pos = pos
        self.refBase = refBase

    def addRead(self, read):
        """
        """

###################################################################################################

import pysam
import sys
import heapq
from collections import defaultdict
from operator import attrgetter
from itertools import tee

###################################################################################################

class Variant(object):
    """
    Simple class to encapsulate variant information.
    """
    def __init__(self, chrom, pos, ref, alt):
        """
        Constructor.
        """
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt

    def __hash__(self):
        """
        So this object can be stored in a set or dictionary
        """
        return hash( (self.chrom, self.pos, self.ref, self.alt) )

    def __eq__(self, other):
        """
        Comparison for equality
        """
        return self.chrom == other.chrom and self.pos == other.pos and self.ref == other.ref and self.alt == other.alt

    def __cmp__(self, other):
        """
        """
        return cmp(self.chrom, other.chrom) or cmp(self.pos, other.pos) or cmp(self.ref, other.ref) or cmp(self.alt, other.alt)

    def __repr__(self):
        """
        For printing.
        """
        return "%s:%s %s --> %s" %(self.chrom,self.pos,self.ref,self.alt)

    def __str__(self):
        """
        For conversion to string.
        """
        return self.__repr__()

###################################################################################################

class Node(object):
    """
    Encapsulate information related to a node of the de-Bruijn graph,
    i.e. a sequence k-mer and related data.
    """
    def __init__(self, sequence, colour, position, weight):
        """
        Construct a Node with a k-mer of sequence and a colour.
        """
        self.sequence = sequence
        self.colours = set([colour])
        self.position = position
        self.weight = weight

    def __hash__(self):
        """
        Allows the Node to be hashed.
        """
        return hash( (self.sequence,) )

    def __cmp__(self, other):
        """
        Allows sorting of nodes by position.
        """
        return cmp(self.position, other.position)

    def __eq__(self, other):
        """
        General comparison function
        """
        return (self.sequence == other.sequence)

    def getPrefix(self):
        """
        Return the prefix sequence.
        """
        return self.sequence[:-1]

    def getSuffix(self):
        """
        Return the suffix sequence.
        """
        return self.sequence[1:]

    def __repr__(self):
        """
        Returns a string representation of the node.
        """
        return "%s (%s). Pos = %s" %(self.sequence, self.colours, self.position)

    def __str__(self):
        """
        Allows node to be converted into a string.
        """
        return self.__repr__()

###################################################################################################

class Edge(object):
    """
    Encapsulates an edge of the de-Bruijn graph.
    """
    def __init__(self, startNode, endNode, colours, weight):
        """
        """
        self.startNode = startNode
        self.endNode = endNode
        self.colours = colours
        self.weight = weight

    def __repr__(self):
        """
        Returns a string representation of the edge.
        """
        return "%s --> %s (colour = %s, weight = %s)" %(self.startNode, self.endNode, self.colours, self.weight)

    def __str__(self):
        """
        Allows edge to be converted into a string.
        """
        return self.__repr__()

###################################################################################################

class DeBruijnGraph(object):
    """
    """
    def __init__(self, kmerSize):
        """
        De-Bruijn graph constructor. Create a load of empty cache variables.
        """
        self.kmerSize = kmerSize
        self.edgesByNode = defaultdict(list)
        self.nodesByPrefix = defaultdict(list)
        self.nodesBySuffix = defaultdict(list)
        self.nodes = {}

    def nodesAreConnected(self, startNode, endNode):
        """
        Return True if a given pair of nodes have an edge between them,
        and False otherwise.
        """
        theEdges = self.edgesByNode[startNode]

        for edge in theEdges:
            print edge
            if edge.endNode == endNode:
                return True

        return False

    def addNode(self, node):
        """
        """
        if node not in self.nodes:
            self.nodes[node] = node

            # Add edges relevant to this node. Any node which has the same prefix
            # as the suffix of this node is connected to it, if they are the same colour.
            for thisNode in self.nodesByPrefix[node.getSuffix()]:
                self.edgesByNode[node].append(Edge(node, thisNode, node.colours, 1))

            for thisNode in self.nodesBySuffix[node.getPrefix()]:
                self.edgesByNode[thisNode].append(Edge(thisNode, node, node.colours, 1))

            self.nodesByPrefix[node.getPrefix()].append(node)
            self.nodesBySuffix[node.getSuffix()].append(node)

        else:
            # Update colours of nodes already in graph.
            self.nodes[node].colours.update(node.colours)
            self.nodes[node].weight += 1

            # Update weight for edges relevant to this node
            for edge in self.edgesByNode[node]:
                edge.weight += 1
                edge.colours.update(node.colours)

###################################################################################################

def breadthFirstSearch(theGraph, sourceNode):
    """
    Search through all nodes in the graph, and build a tree of all nodes which are
    reachable from sourceNode.
    """
    for node in theGraph.nodes:

        if node == sourceNode:
            continue
        elif node.colour != sourceNode.colour:
            continue
        else:
            node.bfsColour = "w"
            node.bfsDistance = None
            node.bfsPredecessor = None

    sourceNode.bfsColour = "g"
    sourceNode.bfsDistance = 0
    sourceNode.bfsPredecessor = None

    queue = [sourceNode]

    while len(queue) > 0:
        node = queue.pop()

        for edge in theGraph.edgesByNode[node]:
            adjacentNode = edge.endNode

            if adjacentNode.bfsColour == "w":
                adjacentNode.bfsColour = "g"
                adjacentNode.bfsDistance = node.bfsDistance + 1
                adjacentNode.bfsPredecessor = node
                queue.append(adjacentNode)

        node.bfsColour = "b"

###################################################################################################

def getAllPathsThroughGraph(theGraph, minWeight):
    """
    Find and return all valid paths through the graph. The total number of
    paths will increase exponentially (?) with connectivity (TODO: check this
    properly). The parameter 'minWeight' specifies the minimum weight required
    for any edge: if weight is below this then the path is not considered.
    """
    posGetter = attrgetter("position")
    sortedNodes = sorted(theGraph.nodes, key=posGetter)
    firstNode = sortedNodes[0]
    lastNode = sortedNodes[-1]
    nodeIndex = 0
    nNodes = len(theGraph.nodes)

    stack = [ [firstNode] ]

    while len(stack) > 0:

        pathSoFar = stack.pop()
        endSoFar = pathSoFar[-1]

        # Got to the end of the path. Add path to list.
        if endSoFar == lastNode:
            #yield pathSoFar[0].position,"".join([x.sequence[0] for x in pathSoFar])
            yield pathSoFar
            startPos = pathSoFar[-1].position
        else:
            for edge in theGraph.edgesByNode[endSoFar]:
                if edge.weight >= minWeight or 0 in edge.colours:
                    nextNode = edge.endNode
                    stack.append( pathSoFar + [nextNode])

###################################################################################################

def depthFirstSearch(theGraph):
    """
    An implementation of depth-first-search taken from Introduction To Algorithms
    by Cormen, Leiserson, Rivest and Stein.

    Here this is implemented as a generator function, yielding successive nodes
    in the depth-first path.
    """
    posGetter = attrgetter("position")
    sortedNodes = sorted(theGraph.nodes, key=posGetter)

    for node in sortedNodes:
        node.dfsColour = "w"

    firstNode = sortedNodes[0]
    theStack = [firstNode]

    while len(theStack) > 0:
        theNode = theStack.pop()

        if theNode.dfsColour == "w":
            yield theNode

            for edge in theGraph.edgesByNode[theNode]:
                theStack.append(edge.endNode)

            theNode.dfsColour = "b"

###################################################################################################

def printAlignments(refStart, refSeq, nodes):
    """
    """
    kmerSize = len(nodes[0].sequence)
    readStart = nodes[0].position + kmerSize//2
    readSeq = "".join([x.sequence[kmerSize//2] for x in nodes])
    readStartInRef = readStart - refStart
    alignmentChars = [" "] * readStartInRef

    for refChar,readChar in zip(refSeq[readStartInRef:], readSeq):

        if refChar == readChar:
            alignmentChars.append(".")
        else:
            alignmentChars.append(readChar)

    print ""
    print ""
    print refSeq
    print "".join(alignmentChars)
    print ""
    print ""

###################################################################################################

#def detectCyclesInGraph(theGraph, minWeight):
#    """
#    Return True if we find any cycles in the graph, and
#    False otherwise. This uses a depth-first search, and
#    works in O(E) using max memory O(V) (I think...).
#    """
#    posGetter = attrgetter("position")
#    sortedNodes = sorted(theGraph.nodes, key=posGetter)
#
#    for node in sortedNodes:
#        node.dfsColour = "w"
#
#    firstNode = sortedNodes[0]
#    theStack = [firstNode]
#
#    while len(theStack) > 0:
#        theNode = theStack.pop()
#
#        if theNode.dfsColour == "w":
#            for edge in theGraph.edgesByNode[theNode]:
#                theStack.append(edge.endNode)
#
#            theNode.dfsColour = "b"
#
#        # Found back-edge, i.e. cycle
#        else:
#            return True
#
#    # DFS completed. No cycles.
#    return False

###################################################################################################

def detectCyclesInGraph(theGraph, minWeight):
    """
    """
    posGetter = attrgetter("position")
    sortedNodes = sorted(theGraph.nodes, key=posGetter)

    for node in sortedNodes:
        node.dfsColour = "w"

    sourceNode = sortedNodes[0]
    destinationNode = sortedNodes[-1]
    theStack = [sourceNode]

    while len(theStack) > 0:
        thisNode = theStack.pop()

        if thisNode == destinationNode:
            return False

        for edge in theGraph.edgesByNode[thisNode]:

            #if edge.dfsColour != "w":
            #    continue

            nextNode = edge.endNode

            if nextNode.dfsColour == "w":
                nextNode.dfsColour = "g"
                theStack.append(nextNode)

            # Found a cycle
            elif nextNode.dfsColour == "g":
                return True

            else:
                pass

        thisNode.dfsColour = "b"

    # Got to the end of the path, without actually reaching the
    # destination node. This will happen if graph is not fully
    # connected.
    return True

###################################################################################################

def buildSequenceFromGraph(theGraph):
    """
    Builds a sequence by walking through the graph following edges in a simple breadth-first
    search manner.
    """
    # Now perform a breadth-first search through the graph.
    thisHap = []
    theseColours = []

    for node in depthFirstSearch(theGraph):

        if node is not None and type(node) != type("-"):
            thisHap.append(node.sequence[0])
            theseColours.append(",".join([str(x) for x in node.colours]))
        else:

            if node is None:
                if len(thisHap) >= kmerSize:
                    yield "".join(thisHap), "".join(theseColours)
                thisHap = []
                theseColours = []

            elif node == "-":
                thisHap.append(node)
            else:
                raise StandardError, "Unexpected value %s" %(node)

###################################################################################################

def getAllKmers(theString, N):
    """
    Return a list of all kmers of size N from the
    given string
    """
    kmers = []
    length = len(theString)

    for x in range(length-N):
        newKmer = theString[x:x+N]
        kmers.append(newKmer)

    return kmers

###################################################################################################

def getAllKmersFromRead(theRead, N):
    """
    Return a list of all kmers of size N from the
    given read.
    """
    theSeq = theRead.seq
    theQuals = [ord(x) - 33 for x in theRead.qual]
    kmers = []
    length = theRead.rlen

    for x in range(length-N):
        newKmer = theSeq[x:x+N]
        newQuals = theQuals[x:x+N]
        kmers.append( (newKmer, min(newQuals) ) )

    return kmers

###################################################################################################

def loadReferenceIntoGraph(theGraph, refFile, chrom, start, end):
    """
    Load k-mers from the specified reference sequence into the
    graph.
    """
    for index,kmer in enumerate(getAllKmers(refSeq, kmerSize)):
        theGraph.addNode(Node(kmer, 0, refStart+index, 1))

###################################################################################################

def loadBAMDataIntoGraph(theGraph, bamFile, chrom, start, end, qualThreshold):
    """
    Load k-mers from the specified BAM file into the graph. K-mers containing
    Ns are ignored, as are k-mers containing low-quality bases.
    """
    for read in bamFile.fetch(chrom, start, end):
        for index,(kmer,minQual) in enumerate(getAllKmersFromRead(read, kmerSize)):
            if minQual >= qualThreshold and kmer.count("N") == 0:
                theGraph.addNode(Node(kmer, 1, read.pos+index, minQual))

###################################################################################################

def getEdgeWeightHistogram(theGraph):
    """
    Return a dictionary of number of edges by weight. This is used to determine
    a threshold for cleaning the graph.
    """
    edgesByWeight = defaultdict(int)

    for node in theGraph.nodes:
        for edge in theGraph.edgesByNode[node]:
            edgesByWeight[edge.weight] += 1

    return edgesByweight

###################################################################################################

def findPathsOnlyInReads(theGraph, minWeight):
    """
    Returns a list of all paths through the graph which include segments only present in
    the read k-mers. This should detect SNPs and Indels, including large deletion breakpoints.
    We can use the position information to spot jumps across the reference in each path.
    """
    bubblePaths = []

    # Find bubbles by looking in best paths
    for path in findNBestPathsThroughGraph(theGraph, 10):
    #for path in getAllPathsThroughGraph(theGraph, minWeight):
        insideLeftFlank = False
        inBubble = False
        bubbleNodes = []
        lastNode = None

        for node in path:

            if not insideLeftFlank:
                if 0 in node.colours and 1 in node.colours:
                    insideLeftFlank = True
                continue

            # Inside flank. Are we still in ref and read?
            elif not inBubble:

                # Ref and read are the same here
                if 0 in node.colours and 1 in node.colours:
                    pass

                # Ref and read are now diverging. Here we are only interested in
                # paths which are in the reads but not the reference.
                if 1 in node.colours and 0 not in node.colours:
                    inBubble = True
                    bubbleNodes.append(lastNode)
                    bubbleNodes.append(node)

            # Currently in bubble
            else:
                # Exiting bubble and returning to ref/read being the same
                if 0 in node.colours and 1 in node.colours:
                    inBubble = False
                    bubbleNodes.append(node)
                    bubblePaths.append(bubbleNodes)
                    bubbleNodes = []

                # Still in bubble. Keep adding to bubbleNodes
                elif 1 in node.colours and 0 not in node.colours:
                    bubbleNodes.append(node)

                # Something odd happened. We went from a path only present in the reads
                # to one only present in the reference. Could be some kind of complex event?
                # Deletion + insertion? Ignore for now.
                else:
                    insideLeftFlank = False # Do nothing until we get back to the reference
                    inBubble = False
                    bubbleNodes = []

            lastNode = node

    return bubblePaths

###################################################################################################

def printGraphInfo(theGraph):
    """
    Print some basic information about the graph. For debug purposes.
    """
    print "There are %s nodes in the graph" %(len(theGraph.nodes))
    print "There are %s edges in the graph" %(sum(len(value) for key,value in theGraph.edgesByNode.iteritems()))

###################################################################################################

def extractVarsFromBubblePaths(theGraph, bubblePath, refSeq, chrom, refStart, refEnd, kmerSize):
    """
    """
    bubbleStartPos = bubblePath[0].position
    bubbleEndPos = bubblePath[-1].position
    bubbleLen = len(bubblePath)
    thisRefSeq = refSeq[bubbleStartPos - refStart: bubbleEndPos - refStart + 1]
    thisReadSeq = "".join([x.sequence[0] for x in bubblePath if 1 in x.colours])

    while len(thisReadSeq) > 0 and len(thisRefSeq) > 0:

        if thisRefSeq[0] == thisReadSeq[0]:
            bubbleStartPos += 1
            thisRefSeq = thisRefSeq[1:]
            thisReadSeq = thisReadSeq[1:]
        else:
            break

    varLen = len(thisReadSeq) - len(thisRefSeq)
    varStartPos = bubbleStartPos

    if varLen == 0:
        thisRefSeq = thisRefSeq[:-1]
        thisReadSeq = thisReadSeq[:-1]

    return Variant(chrom, varStartPos, thisRefSeq, thisReadSeq)

###################################################################################################

def findRegionsWithMissingWeight(theGraph):
    """
    Find and return a list of all contiguous regions in the reference sequence which are not
    covered by the reads.
    """
    posGetter = attrgetter("position")
    sortedNodes = sorted(theGraph.nodes, key=posGetter)

    missingPositions = []
    startRegion = None
    endRegion = None

    for node in sortedNodes:
        if node.weight < 2:

            if startRegion is None:
                startRegion = node.position
                endRegion = node.position
            else:
                endRegion = node.position

        else:
            if startRegion is None:
                continue
            else:

                if endRegion - startRegion > theGraph.kmerSize:
                    missingPositions.append( (startRegion, endRegion) )

                startRegion = None
                endRegion = None

    return missingPositions

###################################################################################################

def extractDeletionsFromReadPaths(refSeq, chrom, refStart, refEnd, kmerSize, minWeight):
    """
    Positions for each node are the start-positions of that k-mer in the reference, if the node
    is either only in the reference or in both reference and reads. So, if we trace all paths
    followed by the read k-mers, we can spot deletions by looking at cases where the position
    jumps by > 1.
    """
    for path in getAllPathsThroughGraph(theGraph, minWeight):

        startPosInRef = None
        startIndexInRef = None

        for index,node in enumerate(path):

            if index > 0:
                diff = node.position - path[index-1].position
                if diff > 1:
                    print diff

            #if startPosInRef is None:
            #    if 1 in node.colours and 0 in node.colours:
            #        startPosInRef = node.position
            #        startIndexInRef = index
            #    continue

            #else:
            #    delEnd = node.position
            #    delStart = path[index-1].position
            #    delLen = delEnd - delStart

            #    if delLen > 1:
            #        print "Found deletion of length %s at %s:%s-%s" %(delLen, chrom, delStart, delEnd)

###################################################################################################

def findNBestPathsThroughGraph(theGraph, N):
    """
    This function uses a modified version of Dijkstra's algorithm for finding the shortest
    path between 2 nodes in a graph. The modified version is written from an outline from
    Wikipedia.
    """
    posGetter = attrgetter("position")
    sortedNodes = sorted(theGraph.nodes, key=posGetter)

    bestPaths = []
    theHeap = []
    countu = defaultdict(int)
    sourceNode = sortedNodes[0]
    destinationNode = sortedNodes[-1]
    heapq.heappush(theHeap, (1, [sourceNode]) )

    while len(theHeap) > 0 and countu[destinationNode] < N:
        Pu = heapq.heappop(theHeap)
        endNodeThisPath = Pu[1][-1]
        countu[endNodeThisPath] += 1

        if endNodeThisPath == destinationNode:
            bestPaths.append(Pu)
            continue

        else:
            if countu[endNodeThisPath] <= N:
                for edge in theGraph.edgesByNode[endNodeThisPath]:
                    Pv = Pu[1] + [edge.endNode]
                    newWeight = 1.0/Pu[0] + edge.weight
                    heapq.heappush(theHeap, (1.0/newWeight, Pv) )

    return [x[1] for x in sorted(bestPaths)]

###################################################################################################

bamFile = pysam.Samfile(sys.argv[1], 'rb')
refFile = pysam.Fastafile(sys.argv[2])

# Test 1kb sized deletion
# Smallish deletion positions on chrom 1
#1:58742911-58745811
#1:105254350-105257189
#1:207291502-207294208

# Test SNP
# A --> G SNP here at 1115529
#chrom = "20"
#start = 1115429
#end = 1115629
#refStart = start - 200
#refEnd = end + 200

# Test small deletion
# Small deletion at 20:1118278, GAA -> G

minBaseQual = 20 # TODO: Make this a function of kmer size? Otherwise SNP calling around repeats is tricky.
windowSize = 1500
minWeight = 2 # I want at least 2 reads supporting any bubble.

#or start in range(1000000, 2000000, 1000):
#for start in range(1118000, 1119000, windowSize//2):
#for start in range(58743000, 58745000, windowSize//2):
#for start in range(207292000, 207294000, windowSize//2):

# Cycles here
#for start in range(80140576, 80142076, windowSize//2):
#chrom = "20"
#for start in range(58743000, 58745000, windowSize//2):
chrom = "1"
theVars = set()

#for start in range(1000000, 1100000, 1000):
#for start in range(1065496, 1065696, 1000):
for line in open('/home/rimmer/Analysis/LargeDeletions/smallishDeletions.txt'):

    chrom = line.split(":")[0]
    start = int(line.split(":")[1].split("-")[0])

    kmerSize = 35
    end = start + windowSize
    refStart = start - 200
    refEnd = end + 200

    print "########################################################################################"
    print "Processing region %s:%s-%s" %(chrom, start, end)
    print "########################################################################################"

    refSeq = refFile.fetch(chrom, refStart, refEnd)

    # Make graph
    theGraph = DeBruijnGraph(kmerSize)
    loadReferenceIntoGraph(theGraph, refFile, chrom, refStart, refEnd)
    loadBAMDataIntoGraph(theGraph, bamFile, chrom, start, end, minBaseQual)

    printGraphInfo(theGraph)

    while detectCyclesInGraph(theGraph, minWeight):
        if kmerSize > 50:
            break
        else:
            print "Found cycles in region %s:%s-%s with kmer size %s. Trying again with kmer size %s" %(chrom, start, end, kmerSize, kmerSize+5)
            kmerSize += 5
            theGraph = DeBruijnGraph(kmerSize)
            loadReferenceIntoGraph(theGraph, refFile, chrom, refStart, refEnd)
            loadBAMDataIntoGraph(theGraph, bamFile, chrom, start, end, minBaseQual)
    else:
        # Find SNPs
        for bubblePath in findPathsOnlyInReads(theGraph, minWeight):
            #print ""
            #print "Found bubble..."
            #for node in bubblePath:
            #    print node
            #print "End of bubble"
            #print ""
            theVars.add(extractVarsFromBubblePaths(theGraph, bubblePath, refSeq, chrom, refStart, refEnd, kmerSize))
            #printAlignments(refStart, refSeq, bubblePath)

        ##extractDeletionsFromReadPaths(refSeq, chrom, refStart, refEnd, kmerSize, minWeight)
        #for rStart,rEnd in findRegionsWithMissingweight(theGraph):
        #    print "Found region %s:%s-%s of size %s with no weight" %(chrom, rStart, rEnd, rEnd - rStart)

for var in sorted(theVars):
    print var




##########################################
#TODO:
#
#3) Alternative deletion caller
#4) Alternative insertion caller
#5) Improvements for small variant calling (smaller k-mer size etc).

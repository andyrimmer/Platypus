"""
Python module for handling intervals. We define an interval as a pair of
values, stored in a tuple (low, high). The interval is half-open, i.e. the
lowest point (low) is contained in the interval, but the highest point (high)
is not.

This module provides functions for storing and manipulating intervals.
"""

###################################################################################################

class Node(object):
    """
    The Node class represents a single node of an interval
    tree. The Node has 5 attributes: 'colour', 'key', 'left', 'right'
    and 'parent'.
    """
    __slots__ = ("colour", "key", "left", "right", "parent")

    def __init__(self, colour, key, left, right, parent):
        """
        Constructor. Sets attributes.
        """
        self.colour = colour
        self.key = key
        self.left = left
        self.right = right
        self.parent = parent

###################################################################################################

class IntervalTree(object):
    """
    The interval tree is an augmented red-black tree. A red-black tree must
    have the following properties:

    1) Every node is either red or black

    2) The root is black

    3) Every leaf is black (and also every leaf has value None)

    4) If a node is red, both its children are black

    5) For each node, all simple paths from the node to descendant roots contain
       the same number of black nodes
    """
    __slots__ = ("root")

    def __init__(self):
        """
        Constructor. Does nothing for now.
        """
        pass

    def leftRotate(self, x):
        """
        Perform a left rotation around the specified node.
        """
        y = x.right
        x.right = y.left

        if y.left is not None:
            y.left.parent = x

        y.parent = x.parent

        if x.parent is None:
            self.root = y
        elif x == x.parent.left:
            x.parent.left = y
        else:
            x.parent.right = y

        y.left = x
        x.parent = y

    def rightRotate(self, x):
        """
        Perform a right rotation around the specified node.
        """
        y = x.left
        x.left = y.right

        if y.right is not None:
            y.right.parent = x

        y.parent = x.parent

        if x.parent is None:
            self.root = y
        elif x == x.parent.left:
            x.parent.left = y
        else:
            x.parent.right = y

        y.left = x
        x.parent = y

    def insert(self, theNode):
        """
        Add a node to the interval tree.
        """
        pass

    def fixup(self):
        """
        Ensure that the tree respects the red-black tree
        properties.
        """
        pass

    def findOverlappingNodes(self, theInterval):
        """
        Returns a sorted list of intervals from all nodes that overlap
        the specified interval.
        """
        pass
###################################################################################################

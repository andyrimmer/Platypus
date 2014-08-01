#!/bin/bash

# Tag a new release of Platypus. This copies (i.e. links) the current svn revision into the tags directory. As this is only a cheap
# copy, it is fast, and uses little space in the repository. The relevant revision can then be checked out directly from tags using
# the command shown below.
#
# svn co file:///home/rimmer/Groups/bioinformatics/svn/repos/tags/${version}

version=$1
echo "Tagging the ${version} release of Platypus"

# I'll keep this commented for now, so I don't create any tags by mistake.
#echo "This is actually commented at the moment"
svn copy file:///home/rimmer/Groups/bioinformatics/svn/repos/trunk file:///home/rimmer/Groups/bioinformatics/svn/repos/tags/${version} -m "Tagging the ${version} release of Platypus"

svn copy file:///home/rimmer/Groups/bioinformatics/svn/repos/trunk/Platypus file:///home/rimmer/Groups/bioinformatics/svn/repos/branches/Platypus2 -m "Branching Platypus"

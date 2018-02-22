from collections import defaultdict
ROOT_END_INDEX = -1
LEAF_END_INDEX = -2
NODE_COUNTER = 0

class SuffixNode:
    def __init__(self, start, end, link):
        global NODE_COUNTER
        self.id = NODE_COUNTER 
        NODE_COUNTER += 1
        self.children = {} # dict of {<edge char>: SuffixNode}
        self.link = link # suffixLink, default None

        #edge data connecting child to parent
        self.start = start # start index. edge stored in the child
        # end index determines the type of node that we're in
        # if root, end index is -1 (end isn't relevant since edge DNE)
        # if leaf, end index is -2, which means the actual end is the endIndex
        # if split, end is the actual split value
        self.end = end
        self.freq = -1

        self.mark_for_removal=False # used to prune out non-maximal suffixes

    def getEnd(self, finalPos):
        e = self.end
        if e == LEAF_END_INDEX:
            e = finalPos
        return e
    
    def edgeLength(self, finalPos):
        """ finalPos is the end of the string thus far"""
        return self.getEnd(finalPos) - self.start + 1
    
    def toStr(self, finalPos, text):
        e = self.getEnd(finalPos)
        s = text[self.start : e+1]
        return "Node: %s, start: %s, end: %s, text: %s, freq: %s" % (self.id, self.start,e, s, self.freq)

    def isLeaf(self):
        return len(self.children) == 0

class SuffixTree:
    def __init__(self, text):
        """
        text is the input string
        """

        # root is the start of the tree. Has no parent
        self.root = SuffixNode(-1, ROOT_END_INDEX, None)

        self.text = text
        # lastNode is the last node that was created
        # we need to keep track of this node since we might
        # need to adjust it's suffix link
        self.lastNode = None

        # below are the triple that represents each new step
        # activeNode is the currently active node according to the algorithm
        # initialized as root
        self.activeNode = self.root
        # activeIndex represents the currently active edge
        self.activeEdge = -1
        # how many levels down from the activeIndex do we care about
        self.activeLength = 0
        # remainder: how many suffixes are pending
        self.remainder = 0

        # endIndex is the SharedIndex for which the corresponding 
        # suffix was the last leaf added
        self.endIndex = -1

    # --------- MAIN INTERFACE ---------------
    def PrintTree(self):
        self.PrintTreeHelper(self.root, 0)
    
    def BuildTree(self):
        # Set up root
        # for each letter in the text, build the suffix tree
        for i in range(len(self.text)):
            self.AdvanceTree(i)
        self.SetFrequencyCounts(self.root)
    
    def PruneTree(self):
        """
        TODO, can do this in the same pass as the SetFrequencyCounts
        For now, keep all motifs that have length > 1 and freq > 1
        """
        self.PruneTreeHelper(self.root,None, None)
    
    # --------- PRIVATE

    def PrintTreeHelper(self, currNode, spaces):
        spacesStr = "\t"*spaces
        if currNode is None:
            return
        print spacesStr + "{"
        
        if currNode.start is not -1: # not a root
            end = currNode.getEnd(self.endIndex)        
            s = self.text[currNode.start : end+1]
            print spacesStr + "\t%s" % currNode.toStr(self.endIndex, self.text)
            # print spacesStr + "\tNode: %s, start: %s, end: %s, text: %s" % (currNode.id, currNode.start, end, s)
        else:
            print "\t root, Node %s" % currNode.id
        for c, n in currNode.children.iteritems():
            self.PrintTreeHelper(n, spaces + 1)
        
        print spacesStr + "}"
    
    def SetFrequencyCounts(self, node):
        if node.freq != -1:
            return node.freq # already computed frequency on this node
        # compute frequency
        if node.isLeaf():
            node.freq = 1
        else:
            f = 0
            for _, c in node.children.iteritems():
                f += self.SetFrequencyCounts(c)
            node.freq = f

        # if there is non-root suffix link from this node, compute the suffix link
        # pointee's frequency tree
        if node.link is not None and node.link is not self.root:
            self.SetFrequencyCounts(node.link)
            if node.link.freq == node.freq: # the pointee to the suffix link can be removed
                node.link.mark_for_removal = True
        node.link = None
        return node.freq

    def PruneNode(self, node, prev, prev_edge):
        node.children = None
        node.link = None
        prev.children.pop(prev_edge)


    def PruneTreeHelper(self, node, prev, prev_edge):
        """
        Prune the subtree governed by this node. Also get rid of all suffix links
        @params node: the node who starts the subtree to be pruned
                prev: the parent of node (None if root)
                e: the edge connecting prev and Node
        """
        if node.mark_for_removal: # non-maximal tagged in previous pass
            self.PruneNode(node, prev, prev_edge)
            return

        if node.freq == 1 and prev is not None:
            # should prune since only has one occurrence
            self.PruneNode(node, prev, prev_edge) # prune from parent
            return
        
        # prune the subtree
        for e, c in node.children.items(): # items over iteritems since deleting from dict
            self.PruneTreeHelper(c, node, e) # try to prune children
        
        if node.edgeLength(self.endIndex) == 1 and prev is self.root: # single character motif from root
            # if the node has a single child, pull the child inward
            if len(node.children) == 1:
                e, c = node.children.items()[0]
                c.start = node.start # make c a more general version of node
                self.PruneNode(node, prev, prev_edge)
                prev.children[prev_edge] = c # replace node with its child
            elif len(node.children) == 0:
                self.PruneNode(node, prev, prev_edge)
        
    def AdvanceTree(self, cursor):
        """
            Advance the tree by moving onto the next index
        """
        # SETUP
        self.endIndex = cursor # update the last leaf
        self.remainder += 1
        self.lastNode = None

        # add the remainder suffixes
        while self.remainder > 0:
            if (self.activeLength == 0):
                self.activeEdge = cursor # want to insert at root
            activeChar = self.text[self.activeEdge] # character at activeEdge

            # Check to see if there is already a subtree
            # with activeEdge from activeNode
            if activeChar not in self.activeNode.children:
                # there is no active edge already at the active node
                # so we create a new leaf edge
                newLeaf = SuffixNode(cursor, LEAF_END_INDEX, self.root)
                self.activeNode.children[activeChar] = newLeaf

                # set the suffix pointer if there is an old node that needs
                # to be set to the current active node
                if self.lastNode is not None:
                    self.lastNode.link = self.activeNode 
                    self.lastNode = None
            else:
                # there exists an appropriate edge coming out of active node
                # perform skip/count optimization
                nextNode = self.activeNode.children[activeChar]
                nextEdgeLength = nextNode.edgeLength(self.endIndex)
                # nextNode.freq[self.activeLength] += 1
                if self.activeLength >= nextEdgeLength:
                    # Skip/Count trick. We can skip to the next node
                    self.activeEdge += nextEdgeLength
                    self.activeLength -= nextEdgeLength
                    self.activeNode = nextNode
                    continue

                currProcessedChar = self.text[nextNode.start + self.activeLength]
    
                # if the char that we're already on is the same as the newest char 
                # that we're supposed to be extending, do not create a new node
                # Showstopper: update suffix link and then return out for next iteration
                if (currProcessedChar == self.text[cursor]):
                    # if we have a pending suffix link, set it to the active node
                    if self.lastNode is not None and self.activeNode is not self.root:
                        self.lastNode.link = self.activeNode 
                        self.lastNode = None
                    self.activeLength += 1
                    return 

                # our activePoint is in the middle of the current edge, otherwise 
                # we would have skip counted. Our active point is not at the edge 
                # of the currently processed string. Thus we have to split
                # creating a new leaf and a new internal node
                internalEnd = nextNode.start + self.activeLength - 1
                internalNode = SuffixNode(nextNode.start, internalEnd, self.root)
                # replace nextNode as the child of activeNode
                self.activeNode.children[activeChar] = internalNode 

                # internalNode now has two children, nextNode and a new leaf node
                # we advance nextNode by activeLength and add it as a child
                nextNode.start += self.activeLength
                internalNode.children[self.text[nextNode.start]] = nextNode
                # add a new leaf node
                newLeafNode = SuffixNode(cursor, LEAF_END_INDEX, self.root)
                internalNode.children[self.text[cursor]] = newLeafNode

                # since we inserted a new internal node, fix suffix from previous 
                # if we haven't already
                if self.lastNode is not None:
                    self.lastNode.link = internalNode
                
                self.lastNode = internalNode

            # added a suffix so decrement remainder
            self.remainder -= 1
            if self.activeNode is self.root and self.activeLength > 0:
                self.activeLength -= 1
                self.activeEdge = cursor - self.remainder + 1
            elif self.activeNode is not self.root:
                self.activeNode = self.activeNode.link

# abcabxabcd$
s = SuffixTree("abcabxabcd$")
s.BuildTree()
s.PrintTree()
print "----------"
s.PruneTree()
s.PrintTree()
    



import utility as util

class gnode:
    in_edge = []
    out_edge = []

    def __init__(self, sequence):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def getOverlap(self, gnode):
        return util.overlap(self.sequence, gnode.sequence, min_length=3)

    def getBalance(self):
        return len(out_edge) - len(in_edge)

    def isBalanced(self):
        return len(in_edge) == len(out_edge)

    def isSemiBalanced(self):
        return abs(len(in_edge) - len(out_edge)) < 2

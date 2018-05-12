import utility as util

class gnode:
    def __init__(self, sequence):
        self.seq = sequence
        self.in_edge = []
        self.out_edge = []

    def __len__(self):
        return len(self.seq)

    def getOverlap(self, gnode):
        return util.overlap(self.seq, gnode.sequence, min_length=3)

    def getBalance(self):
        return len(out_edge) - len(in_edge)

    def isBalanced(self):
        return len(in_edge) == len(out_edge)

    def isSemiBalanced(self):
        return abs(len(in_edge) - len(out_edge)) < 2

    """
    #Burrows-Wheeler transform generation implementation from Wikipedia
    """
    def bwt(self):
        s = self.seq
        """Apply Burrows-Wheeler transform to input string."""
        assert "\002" not in s and "\003" not in s, "Input string cannot contain STX and ETX characters"
        s = "\002" + s + "\003"  # Add start and end of text marker
        table = sorted(s[i:] + s[:i] for i in range(len(s)))  # Table of rotations of string
        last_column = [row[-1:] for row in table]  # Last characters of each row
        return "".join(last_column)  # Convert list of characters into string

    def ibwt(self, r):
        """Apply inverse Burrows-Wheeler transform."""
        table = [""] * len(r)  # Make empty table
        for i in range(len(r)):
            table = sorted(r[i] + table[i] for i in range(len(r)))  # Add a column of r
        s = [row for row in table if row.endswith("\003")][0]  # Find the correct row (ending in ETX)
        return s.rstrip("\003").strip("\002")  # Get rid of start and end markers

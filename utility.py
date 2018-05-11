import os

notification = "well.wav"

# Determine how much the end of 'a' overlaps with the beginning of 'b'
def overlap(a, b, min_length = 3):
    start = 0
    while 1:
        start = a.find(b[:min_length], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

def play(path):
    os.system("start "+path)

def notify():
    os.system("start " + notification)

def fprint(value):
    print("{:,}".format(int(value)))

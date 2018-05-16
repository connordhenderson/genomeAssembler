import os, random
def play(path="well.wav"):
    os.system("start "+path)

def notify():
    os.system("start " + notification)

def fprint(value):
    print("{:,}".format(int(value)))

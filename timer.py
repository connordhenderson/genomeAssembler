import time, sys

class timer:
    def __init__(self,name):
        self.name = name
        self.start = time.time()

    def start(self):
        self.start = time.time()

    def end(self):
        t = time.time() - self.start
        self.start = time.time()
        return t

    def go(self, text=""):
        return time.time() - self.start

    def pend(self):
        t = time.time() - self.start
        print("%s [%.2f]" % (self.name, time.time() - self.start))
        sys.stdout.flush()
        self.start = time.time()
        return t

    def pgo(self, text=""):
        if text != "":
            text = "("+text+")"
        print("%s [%.2f] %s" % (self.name, time.time() - self.start, text))
        sys.stdout.flush()
        return time.time() - self.start

import sys

path = None

for arg in sys.argv:
    if arg.startswith("f="):
        path = arg[2:]

if path != None:
    sequence = ""
    with open(path) as file:
        try:
            origin_line = None
            origin_started = False
            while 1:
                line = next(file)
                if origin_started:
                    origin_line = line.strip().split(" ")[1:]
                    print(origin_line)
                    sequence += "".join(origin_line)
                if (line.startswith("ORIGIN")):
                    origin_started = True
                    print("found origin")
        except StopIteration:
            print(len(sequence))
            print("parsed sequence")
        file.close()

    with open("Data/gbff.dat",'w') as file:
        file.write(sequence.upper())

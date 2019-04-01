from sys import argv


values = None
try:
   values = [ int(i) for i in argv[1:] ]
except:
   exit(1)

print sum(values)


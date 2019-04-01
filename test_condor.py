#!/usr/bin/env python
from sys import argv
import time
import math
values = None
try:
   values = [ int(i) for i in argv[1:] ]
except:
   exit(1)



for i in arange(1.e5):
    x = math.sqrt(i)

time.sleep(100)




print sum(values)


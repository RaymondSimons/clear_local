#!/usr/bin/env python
from sys import argv
import time
values = None
try:
   values = [ int(i) for i in argv[1:] ]
except:
   exit(1)


time.sleep(40)


print sum(values)


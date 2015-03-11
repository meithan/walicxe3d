# Creates a visit master file for the specified outputs

# If executed without arguments, runs in interactive mode.
# Otherwise, the user must provide the number of processors, followed by one
# or several output numbers. An output number can be:
#  1) a single output number, X
#  2) a range of the form X-Y, going from X to Y (inclusive)
#  3) a range of the form X-Y-S, going from X to Y in steps of S

import sys

# ===============

def err_msg():
  print
  print "Input error!"
  print
  print "Program usage:"
  print "master (no arguments: interactive mode)"
  print "    -OR-"
  print "master <processors> <outputspec> [,<outputspec>]"
  print "where <outputspec> is an output number specification:"
  print "1) a single output number, X"
  print "2) a range of the form X-Y, from X to Y (inclusive)"
  print "3) a range of the form X-Y-S, from X to Y in steps of S"
  print "If more than one output number specification is given, they"
  print "must be comma-separated."
  sys.exit()

# ===============

if len(sys.argv) == 1:
  ans = raw_input("Number of processors: ")
  try: nprocs = int(ans)
  except: err_msg()
  outstr = raw_input("Enter output numbers: ")
else:
  try:
    nprocs = int(sys.argv[1])
    outstr = " ".join(sys.argv[2:])
  except:
    err_msg()

# Parse outstr to build list of output numbers
outlist = []
outspecs = map(lambda x: x.strip(), outstr.split(","))
for spec in outspecs:
  if spec.count("-") == 0:
    try: outlist.append(int(spec))
    except: err_msg()
  elif spec.count("-") in [1,2]:
    ns = spec.split("-")
    try:
      a = int(ns[0])
      b = int(ns[1])
      if len(ns) == 3:
        s = int(ns[2])
      else:
        s = 1
      for n in range(a,b+1,s):
        outlist.append(n)
    except:
       err_msg()
  else:
     err_msg()
outlist.sort()
     
# Write master
f = open("master.visit","w")
f.write("!NBLOCKS%6i" % nprocs)
for outnum in outlist: 
  for proc in range(nprocs):
    f.write("\nBlocks%03i.%04i.vtk" % (proc, outnum))
f.close()
print "Wrote master.visit"


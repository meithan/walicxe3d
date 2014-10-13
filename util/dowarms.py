import os

# ==================================

def modify_source (fname, mod_list):

  # Backup file
  os.system("cp %s %s.bak" % (fname,fname))

  # Load file
  f = open(fname, 'r')
  lines = f.readlines()    
  f.close()

  # Modify source
  ncount = 0
  for modif in mod_list:
    line_num = modif[0]
    new_text = modif[1]
    lines[line_num-1] = new_text+"\n"
    ncount += 1

  # Write back modified file
  f = open(fname, 'w')
  f.writelines(lines)
  f.close()

  print "Peformed %i changes in %s" % (ncount,fname)

# ==================================

# Code for automated warm starts

# CONFIG
nmin = 20
nmax = 20
run = False

for nout in range(nmin,nmax+1):

  # ==================

  # Parameters file
  fname = "./source/parameters.f90"
  mod_list = []

  str = "  logical, parameter :: dowarm = .true."
  mod_list.append([47,str])

  str = '  character(*), parameter :: warm_file = "/home/claudio/G352/S10/BIN/State.%04i.dat"' % (nout)
  mod_list.append([50,str])

#  print mod_list
  modify_source(fname, mod_list)

  # ==================

  # Makefile
  fname = "./Makefile"
  mod_list = []

  str = "PROGRAM= G352-S10-%i" % (nout)
  mod_list.append([6,str])

#  print mod_list
  modify_source (fname, mod_list)

  # ==================

  # PBS launcher
  fname ="./launch.pbs"
  mod_list = []
  
  str = "#PBS -N G352-S10-%i" % (nout)
  mod_list.append([3,str])

  str = "#PBS -e G352-S10-%i.err" % (nout)
  mod_list.append([4,str])

  str = "#PBS -o G352-S10-%i.out" % (nout)
  mod_list.append([5,str])

  str = 'EXE="G352-S10-%i"' % (nout)
  mod_list.append([23,str])

#  print mod_list
  modify_source (fname, mod_list)

  # ==================

  # Compile code
  os.system("make clean")
  os.system("make")
    
  # Run code - if asked
  if (run): os.system("qsub launch.pbs")

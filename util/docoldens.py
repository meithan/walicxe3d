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

# Code for automated coldens maps

launch = True

# ==================================

for accel in ["PAR"]:
  for angle in [0]:

    program = "coldens%s%.2i" % (accel,angle)

    # ==================
    # 1) Modify coldens.f90   

    fname = "./source/coldens.f90"
    mod_list = []
 
    str = "  real, parameter :: rot_beta = %3.1f" % (angle)
    line = 93
    mod_list.append([line,str])

    str = "  integer, parameter :: accel = ACCEL_%s" % (accel)
    line = 165
    mod_list.append([line,str])

#    print mod_list
    modify_source(fname, mod_list)

    # ==================
    # 2) Modify Makefile   

    fname = "./Makefile"
    mod_list = [] 

    str = "PROGRAM= %s" % (program)
    line = 6
    mod_list.append([line,str])

#    print mod_list
    modify_source(fname, mod_list)

    # ==================
    # 3) Compile code

    os.system("make clean")
    os.system("make coldens")
    os.system("mv coldens %s" % (program))

    # ==================
    # 4) Modify PBS launch script

    fname = "./launch.pbs"
    mod_list = []
  
    str = "#PBS -N %s" % (program)
    line = 3
    mod_list.append([line,str])

    str = "#PBS -e %s.err" % (program)
    line = 4
    mod_list.append([line,str])

    str = "#PBS -o %s.out" % (program)
    line = 5
    mod_list.append([line,str])

    str = 'EXE="%s"' % (program)
    line = 23
    mod_list.append([line,str])

    print mod_list
    modify_source (fname, mod_list)

    # ==================
    
    # 5) Launch program
    if (launch): os.system("qsub launch.pbs")

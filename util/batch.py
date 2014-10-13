import os

# Source to modify
fname = "./source/coldens.f90"

# Backup file
os.system("cp %s %s.bak" % (fname, fname))

for accel in ["PAR","ISO"]:
  for angle in [0,30,45,50,55,60,65,70,75,90]:

    # Replace accel method
    new_line = "  integer, parameter :: accel = ACCEL_%s" % accel
    line_num = 165
    cmd = "sed '%is/.*/%s/' %s > foo" % (line_num, new_line, fname)
    print cmd
    os.system(cmd)
    
    # Replace angle        
    new_line = "  real, parameter :: rot_beta = %02i.0" % angle
    line_num = 93
    cmd = "sed '%is/.*/%s/' foo > bar" % (line_num, new_line)
    print cmd
    os.system(cmd)

    # Overwrite original file
    cmd = "cp bar source/coldens.f90"
    print cmd
    os.system(cmd)

    # Compile code
    os.system("make coldens")
    
    # Run code
    os.system("./coldens")

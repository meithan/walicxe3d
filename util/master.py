# Creates a visit master file for the specified outputs

# If executed without arguments, runs in interactive mode.
# Otherwise, the user must provide command-line arguments. See --help for details.

import argparse
import sys

# ===============

parser = argparse.ArgumentParser(
  prog="python master.py",
  description="""Creates a visit master file for the specified outputs.

Call without arguments for interactive mode
-OR-
Call with the following command-line arguments:""",
formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument("num_procs", type=int, nargs="?", help="The number of processes")
parser.add_argument("output_spec", nargs="?", help="The output or range of outputs, dash-separated, which must be one of:\n1) A single output number, X\n2) A range of the form X-Y, meaning from X to Y (inclusive)\n3) A range of the form X-Y-S, meaning from X to Y in steps of S")
parser.add_argument("-i", "--interactive", action="store_true", help="Run in interactive mode (default if no arguments passed)")

args = parser.parse_args()

if len(sys.argv) == 1 or args.interactive:
  # Interactive mode
  print("(Running in interactive mode, see --help for other options)")
  nprocs = int(input("Number of processes: "))
  output_spec = input("Output numbers (X, X-Y or X-Y-S): ")
else:
  # Parse command-line arguments
  nprocs = args.num_procs
  output_spec = args.output_spec

if output_spec.count("-") == 0:
  outlist = [int(output_spec)]
elif output_spec.count("-") in [1, 2]:
  tokens = output_spec.split("-")
  start = int(tokens[0])
  end = int(tokens[1])
  if len(tokens) == 3:
    step = int(tokens[2])
    outlist = list(range(start, end+1, step))
  else:
    outlist = list(range(start, end+1))

# Write master file
f = open("master.visit", "w")
f.write("!NBLOCKS {}\n".format(nprocs))
for outnum in outlist: 
  for proc in range(nprocs):
    f.write("Blocks{:03d}.{:04d}.vtk\n".format(proc, outnum))
f.close()
print("Wrote master.visit")


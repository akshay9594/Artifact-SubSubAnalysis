

import os,sys

from ExecScripts import Exp1_Script,Exp2_Script

print("\n*************************************************************************\n")
print("\t\t>Artifact execution for reproducing the results<")


root = os.getcwd()

print("This experiment measures performance improvement for 3 benchmarks")
Exp2_Script.RunExp(root)

print("\n*************************************************************************\n")

sys.exit()

#Main code
#Step 1: Select the benchmark to be evaluated
val = input("\t1. Experiment 1\n\t2. Experiment 2\n\t3. Both Experiments\n->Select the Experiment you want to run by entering the number:\n")


if(val == '1'):
    print("\nRunning Experiment 1...")
    print("This experiment measures performance improvement for 3 benchmarks")
    Exp1_Script.RunExp(root)

elif(val == '2'):
    print("\nRunning Experiment 2...")

elif(val == '3'):
    print("\nRunning both Experiments...")

else:
    print("Invalid Selection")
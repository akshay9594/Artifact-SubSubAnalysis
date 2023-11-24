

import os,sys

import subprocess

from ExecScripts import Exp1_Script,Exp2_Script

from utils import runcmd



print("\n*************************************************************************\n")
print("\t\t>Artifact execution for reproducing the results<")
print("\n*************************************************************************\n")


root = os.getcwd()

input_matrices_direc = root + "/input_matrices/"

if(os.path.exists(input_matrices_direc) == False):
    os.mkdir(input_matrices_direc)

os.chdir(input_matrices_direc)

input_matrix_dict = {'gsm_106857':'https://suitesparse-collection-website.herokuapp.com/MM/Dziekonski/',
               'dielFilterV2clx':'https://suitesparse-collection-website.herokuapp.com/MM/Dziekonski/',
               'af_shell1':'https://suitesparse-collection-website.herokuapp.com/MM/Schenk_AFE/',
               'inline_1':'https://suitesparse-collection-website.herokuapp.com/MM/GHS_psdef/',
               'crankseg_1':'https://suitesparse-collection-website.herokuapp.com/MM/GHS_psdef/',
               'spal_004':'https://suitesparse-collection-website.herokuapp.com/MM/Mittelmann/'}
  

# #Download the required external input matrices.
# #Tar files of the matrices are download. The matrix is extracted and
# # then the tar file is deleted.

print("Downloading the required external input matrices from the SuiteSparse Matrix Collection...")

for matrix in input_matrix_dict.keys():
    download_url = input_matrix_dict[matrix]
    matrix_tar = matrix + '.tar.gz'
   
    if(os.path.exists(input_matrices_direc + matrix) == False):

        download_url = download_url + matrix_tar

        runcmd("wget "+download_url, verbose = True)

        runcmd("tar -xf "+matrix_tar, verbose = True)

        runcmd("rm -d "+matrix_tar, verbose = True)
    
print("\nMatrices downloaded and placed within the input_matrices directory!!")
os.chdir(root)

#Main code
#Step 1: Select the benchmark to be evaluated
# val = input("\t1. Experiment 1\n\t2. Experiment 2\n\t3. Both Experiments\n->Select the Experiment you want to run by entering the number:\n")

print("\nRunning Experiment 1...\n")
Exp1_Script.RunExp(root)

# if(val == '1'):
#     print("\nRunning Experiment 1...")
#     print("==>This experiment measures performance improvement for 3 benchmarks")
#     Exp1_Script.RunExp(root)

# elif(val == '2'):
#     print("\nRunning Experiment 2...")
#     Exp2_Script.RunExp(root)

# elif(val == '3'):
#     print("\nRunning both Experiments...")

# else:
#     print("Invalid Selection")
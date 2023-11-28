

import os

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

print("(a)  Step 1:  Downloading the required inputs ...")

for matrix in input_matrix_dict.keys():
    download_url = input_matrix_dict[matrix]
    matrix_tar = matrix + '.tar.gz'
   
    if(os.path.exists(input_matrices_direc + matrix) == False):

        download_url = download_url + matrix_tar

        runcmd("wget "+download_url, verbose = True)

        runcmd("tar -xf "+matrix_tar, verbose = True)

        runcmd("rm -d "+matrix_tar, verbose = True)
    
print("\n\t-->Downloads Complete!!")
os.chdir(root)

print("\n--------------------------------------------------------------------------\n")

print("(b)  Step 2:  Experiment selection and execution\n")

#Main code
# Step 2: Select the Experiment to be executed
val = input("\t1. Experiment 1 (Measures performance improvement for 3 benchmarks)\n\t" +
            "2. Experiment 2 (Compares the performance impact of the Base technique and New technique)\n\t" +
            "3. Both Experiments\n\n\t-->Select the Experiment you want to run by entering the number: ")

# Run the experiment
if(val == '1'):
    print("\n\nRunning Experiment 1...")
    Exp1_Script.RunExp(root)
    print("\nExperiment finished! Results written to the Reports directory and Graphs to the Graphs directory!!")

elif(val == '2'):
    print("\n\nRunning Experiment 2...")
    Exp2_Script.RunExp(root)
    print("\nExperiment finished! Results written to the Reports directory and Graphs to the Graphs directory!!")

elif(val == '3'):
    print("\n\nRunning Experiment 1...")
    Exp1_Script.RunExp(root)

    print("\n--------------------------------------------------------------------------\n")

    print("\nRunning Experiment 2...")
    Exp2_Script.RunExp(root)
    print("\nExperiment finished! Results written to the Reports directory and Graphs to the Graphs directory!!")

else:
    print("\nInvalid Selection")


print("\n*************************************************************************\n")
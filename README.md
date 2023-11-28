# Artifact Description for the Evaluation of Subscripted Subscript Analysis

[![DOI](https://zenodo.org/badge/670291996.svg)](https://zenodo.org/doi/10.5281/zenodo.10213905)

This README describes how to evaluate the artifact for the paper: 
"Recurrence Analysis for Automatic Parallelization of Subscripted Subscripts" submitted to
the ACM SIGPLAN Symposium on Principles and Practice of Parallel Programming (PPoPP) 2024.

## What is reproduced?
The Artifact reproduces major parts of the evaluation results of Experiment 1 and 
Experiment 2 mentioned in Section 4 of the paper. More specifically, the following results are
reproduced:

### For Experiment 1:
The speedup graphs for benchmarks AMGmk, SDDMM and UA(transf) shown in Figure 14. The scripts
produce one graph for each benchmark. The graphs plot the performance improvement of the
Cetus parallel codes (OpenMP parallelization) v/s the Serial baseline.

Note: 
- In the interest of time, for UA-NAS benchmark, CLASSES A, B, C are used.
- The results in Figure 13 of the paper are NOT reproduced, but the source files are provided
  within the *Cetus-Output-WithoutSubSub* subdirectories.
- The scripts measure and plot the performance improvement for the maximum number of cores 
  available on the machine. The scripts cannot vary the number of available cores.

### For Experiment 2

The speedup graph shown in Figure 17 of the paper. The graph compares the impact of two techniques -
*Cetus+BaseAlgo* and *Cetus+NewAlgo* on the performance of 12 benchmarks listed in Table 1 of the
paper. The table also shows the inputs used for each benchmark. We used MATRIX2, dielFilterV2clx and 
CLASS A as input datasets for the AMGmk, SDDMM and UA applications in this experiment.

Note:
  The scripts measure and plot the performance improvement for the maximum number of cores 
  available on the machine. The scripts cannot vary the number of available cores.

## Prerequisities
### Software
 - Linux (OS tested with : CentOS v7.4, Ubuntu v22.04)
 - GNU C Compiler (GCC) v4.8.5 and above
 - Python v3.8.0 and above
 - OpenMP v4.0 and above
 - awk
 - bc (linux calculator)
 - gfortran

### Python packages required
- Non built-in packages:
1. matplotlib
 - Built-in packages:
2. subprocess
3. re
4. math
5. os
6. Numpy
7. shutil

### Hardware
 - Machine with x86-64 processors (preferably Sky Lake and beyond)
 - ~2.5GB of disk space

## Obtaining the Codes
1. The codes can be obtained from the Zenodo repository using the DOI -
   [![DOI](https://zenodo.org/badge/670291996.svg)](https://zenodo.org/doi/10.5281/zenodo.10213905)
   
   A zipped folder is downloaded.Unzip the folder and place it in the root directory of your system.

2. A docker image of the artifact is also available. 

   - The image is hosted at : 
      
      https://hub.docker.com/r/akshaybhosale9594/artifact-subsubanalysis/tags

   - First start the docker engine on your system and then pull and run the
     docker image using the following commands:

   ```
   docker pull akshaybhosale9594/artifact-subsubanalysis
   ```
   ```
   docker run -t -i akshaybhosale9594/artifact-subsubanalysis:v1.0
   ```
   - Use the *master.py* script to execute the experiments using the instructions below.

## Code description
- The source code files for each experiment are placed in the directories -- *Experiment_1*
    and *Experiment_2*.
- For each benchmark, the Baseline and Optimized (Technique(s)_Applied) source files  
  are provided.
- The optimized files refer to the Cetus translated versions of the original source code.
- In *Experiment_1*, the source files are further arranged according to the inputs to a 
  benchmark if the benchmark uses internally generated (within the code) inputs. 
  E.g. for amgmk, the source files are arranged into directories *MATRIX1* through *MATRIX5* 
  as these matrices are internally generated. 
- In *Experiment_2*, the optimized source files are further arranged into the directories
  *Base_Technique* (referring to the Base Algorithm of [5]) and *New_Technique* 
  (referring to the New Algorithm presented in the paper).

## Installing the dependencies

### Installing awk, bc and gfortran

Following commands can be used to install awk, bc and gfortran:

```
 sudo apt-get install -y gawk
 sudo apt-get install -y bc
 sudo apt-get install -y gfortran
```

### Installing the non built-in python packages
The non built-in packages have been mentioned in requirements.txt and can be installed 
using the python package manager : pip or pip3, using the command:
```
pip install -r requirements.txt
```
Note: The dependencies are pre-installed in the docker container.

## Compiling and Running the Codes

Expected Completion time: ~1 hour for each experiment, total ~2 hrs.
Each benchmark is run 3 times.

### The master script:

- The master script to perform the evaluation is the python script : "master.py".
- The master script is an interactive script and uses user input to determine which 
  experiment to run -- *Experiment 1*, *Experiment 2* or both.
- The script automatically downloads external inputs required to run some benchmarks.
- Run the master script using the following command:

    ```
    python3 master.py
    ```
  If using the docker image, use the command:

    ```
    python master.py
    ```
- The script compiles the benchmark files, executes the binary and records the measurements. 

### Generated Results:

- The master script generates the following:

1. Execution time and Speedup Reports for each benchmark and each experiment
   - The result reports for each experiment are placed in the *Reports* directory.
   - The reports show the average execution times of the baseline and optimized codes
     and the calculated average speedups.

2. Graphs for each Experiment
   - The generated graphs for each experiment are placed in the *Graphs* directory,
     within the *Exp-1* and *Exp-2* subdirectories.
   - *Exp-1* stores the graphs of Experiment 1 and *Exp-2* stores the graphs of Experiment 2.

## Visualizing Graphs in Docker
1. The most convenient way to access and visualize graphs from running
    docker containers is to use the Docker extension within a code editor
    such as Visual Studio Code (VSCode).
2.  After installing the Docker extension in VSCode, use the Remote Ex-
    plorer tab to establish a remote connection with the running container.
3.  The artifact files will be listed under *files* > *root* > *artifact-subsubAnalysis*
     
## Translating an input code through the Cetus executable (Optional):

- A Cetus executable has been provided which can be used to perform sanity checks, ensuring that
  Cetus with subscripted subscript analysis enabled, generates the expected optimized code.
- The Cetus executable can be found in the directory *Cetus-bin*.
- Note that only Serial Baseline source codes can be translated.
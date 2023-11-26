# Artifact Description for the Evaluation of Subscripted Subscript Analysis

This README describes how to evaluate the artifact for the paper: 
"Recurrence Analysis for Automatic Parallelization of Subscripted Subscripts" submitted to
the ACM SIGPLAN Symposium on Principles and Practice of Parallel Programming (PPoPP) 2024.


## What is reproduced?
The Artifact reproduces major parts of the evaluation results of Experiment 1 and 
Experiment 2 mentioned in the paper. More specifically, the following results are
reproduced:

### For Experiment 1
The speedup graphs for benchmarks AMGmk, SDDMM and UA(transf) shown in Figure 14. The scripts
produce one graph for each benchmark. The graphs plot the performance improvement of the
Cetus parallel codes (OpenMP parallelization) v/s the Serial baseline.

Note: The scripts measure and plot the performance improvement for the maximum number of cores 
available on the machine. The scripts cannot vary the number of cores available cores.

### For Experiment 2

The speedup graph shown in Figure 17 of the paper. The graph compares the impact of two techniques -
*Cetus+BaseAlgo* and *Cetus+NewAlgo* on the performance of 12 benchmarks listed in Table 1 of the
paper. The table also shows the inputs used for each benchmark. We used MATRIX2, dielFilterV2clx and 
CLASS A as input datasets for the AMGmk, SDDMM and UA applications in this experiment.

Note: The scripts measure and plot the performance improvement for the maximum number of cores 
available on the machine.

## Prerequisities
### Software
 - Linux (OS tested with : CentOS v7.4, Ubuntu v21.04)
 - GNU C Compiler (GCC) v4.8.5 and above
 - Python v3.8.0 and above
 - OpenMP v4.0 and above

### Python Packages required
- Non built-in packages:
1. matplotlib
- Built-in packages:
2. subprocess
3. re
4. math
5. os
6. Numpy

### Hardware
 - Machine with x86-64 processors (Sky Lake and beyond)

## Obtaining the Codes
The codes can be obtained from the Zenodo repository using the DOI.

## Code Description
- The source code files for each experiment are placed in the directories -- *Experiment 1*
    and *Experiment 2*.
- For each benchmark, the Baseline and optimized (Technique(s)_Applied) source files  
  are arranged.
- The optimized files refer to the Cetus translated versions of the original source code.
- For *Experiment 1*, the source files are further arranged according to the inputs to a 
  benchmark if, the benchmark uses internally generated (within the code) inputs. 
  E.g. for amgmk, the source files are arranged into directories *MATRIX1* through *MATRIX5* 
  as these matrices are internally generated. 
- For *Experiment 2*, the optimized source files are further arranged into the directories
  *Base_Technique* (referring to the Base algorithm of [5]) and *New_Technique* 
  (referring to the New algorithm presented in the paper).


### Installing the Non built-in python packages
The non built-in packages can be installed using the python package manager : pip3

Using the command:
```
pip3 install -r requirements.txt
```

## Compiling and Running the Codes

Expected Completion time: ~1 hour for each experiment depending on the machine

### The master script:

- The master script to perform the evaluation is the python script : master.py in the
  root directory.
- The master script is an interactive script and uses user input to determine which 
  experiment to run -- *Experiment 1*, *Experiment 2* or both.
- The script automatically download external inputs required for some benchmarks.
- Run the master script using the following command:

    ```
    python3 master.py
    ```


### Generated Results:

- The master script generates the following:

1. Execution time and Speedup Reports for each benchmark and each experiment
   - The reports are placed in the *Reports* directory for each experiment.
   - The reports show the execution times of the baseline and optimized codes
     and the speedups.

2. Graphs for each Experiment
   - The generated graphs for each experiment are placed in the *Graphs* directory.
   - For Experiments 1 and 2, the graphs are placed in the *Exp-1* and *Exp-2* 
     subdirectories respectively.


## Translating an input code through the Cetus executable (Optional):

- A Cetus executable has been provided which can be used to perform sanity checks, ensuring that
  Cetus with subscripted subscript analysis enabled, generates the expected optimized code.
- The Cetus executable can be found in the directory *Cetus-bin*.
- Note that only Serial Baseline source codes can be translated.
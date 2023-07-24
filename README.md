# Artifact Description for the Evaluation of Subscripted Subscript Analysis

We have developed a new analysis technique for the automatic parallelization of subsripted
subscript loops. The technique analyzes loops that define and/or modify the subscript array
and determines array properties, which is sufficient to parallelize a class of subscripted
subscripts. This repository mentions benchmark source codes used to evaluate two capabilites
of the technique - (a) determining monotonicity of arrays that store intermittent sequences
and (b) determining monotonicity of multi-dimensional subscript arrays. The technique has
been described in detail in our paper - "Recurrence Analysis for Automatic Parallelization
of Subscripted Subscripts" submitted to the Supercomputing 2023 conference. We present the
experimental setup and steps involved to reproduce the results mentioned in the paper.

## Execution environment
The results described in the above submitted paper have been gathered using the following 
execution environment- a compute node with a 20-core Intel Xeon Gold 6230 processors in 
a dual socket configuration, with a processor base frequency of 2.1 GHz, 27.5MB cache and we 
used upto 8GB of DDR4 memory. The application codes were compiled using GCC v4.8.5 with the 
-O3 optimization flag enabled on CentOS v7.4.1708 and we report the mean of 5 application runs. 
We observed an average run-to-run variation of 1.45% and we used one thread per core.

## Dependencies
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

### Hardware
 - Machine with Intel processors (Sky Lake and beyond)

## Obtaining the Codes
The codes can be obtained from Zenodo from the above mentioned DOI. (or)

The codes can be obtained by cloning the git repository using the commands-
```
cd $HOME or cd ~
git clone (https://github.com/akshay9594/Artifact-Supercomputing-2023.git)
cd Artifact-Supercomputing-2023
```
This repository contains the benchmark codes listed below:

| Code  | Source | Original Source link | 
| ------------- | ------------- | ------------- |
| Amgmk-v1.0  | CORAL Benchmark Codes | (https://asc.llnl.gov/coral-benchmarks)
| UA-NAS | NAS Parallel Benchmarks | (https://github.com/akshay9594/SNU_NPB-1.0.3)


### Installing the Non built-in python packages
The non built-in packages can be installed using the python package manager : pip3

Using the command:
```
pip3 install -r requirements.txt
```

## Directory Tree and Code Descriptions
The directory layout of this Artifact is as follows-

```bash
├── Artifact-Supercomputing-2023
│   ├── amgmk-v1.0
│   │   ├── Baselines
│   │   │   ├── Cetus-Output-WithoutSubSub
│   │   │   ├── Serial
│   │   └── Technique_Applied
│   ├── UA-NAS
│   │   ├── Baselines
│   │   │   ├── Cetus-Output-WithoutSubSub
│   │   │   ├── Serial
│   │   └── Technique_Applied
└── README.md
└── requirements.txt
└── run-exp.py
```
### Baseline and Optimized Codes
- The experiments have two Baselines against which performance improvement is measured:
 1. Serial Code
 2. Cetus Parallelized code (without subscripted subscript analysis enabled)

- The optimized codes are the Cetus parallelized codes (With technique enabled)

***The serial baseline codes have been placed in the "Serial" directory
and the Cetus parallelized codes (without subscripted subscript analysis applied) have 
been placed in the "Cetus-Output-WithoutSubSub" directory for each benchmark.***

***The Cetus parallelized codes with subscripted subscript analysis applied (Optimized 
codes) have been placed in the "Technique_Applied" directory for each benchmark.***

## Compiling and Running the Codes

### The master script:

- The master script to perform the evaluation is the python script : run-exp.py
- The master script is an interactive script and uses user input to:
  1. Determine the benchmark to evaluate.
  2. Set the baseline for the experiment (Serial or Cetus parallelized code 
     without technique enabled).
  3. Determine if the experiment should be run using one or all of the input 
     matrices/classes.
  4. Determine the number of runs (of the baseline code and the optimized application code) 
     to determine the average execution time.
  5. If the user chooses to run the experiment for all of the available inputs, 
     the master script generates a graph showing the performance improvements, 
     similar to the graphs presented in the Evaluation section.

Notes: 
1.    The master script does not vary the number of cores for the experiment as shown in the
      results section of the Supercomputing paper. It uses the maximum available cores. The
      user will have to first set the number of cores to use before running the master
      script.


### Generated Results:

- The master script generates the following results:
    - For the amgmk benchmark:
        1. The Average baseline 'application' execution time
        2. Average 'application' execution time of Cetus Parallelized Code (with technique applied)
        3. Average baseline 'kernel' execution time
        4. Average 'kernel' execution time of Cetus Parallelized Code (with technique applied)
        5. The standard deviation for each of the above mentioned timing results
        6. The performance improvement (Time taken by optimized code/Time taken by the baseline)
        7. Displays the number of threads used

    - For the UA benchmark:
        1. Average baseline execution time of the _transf_() routine
        2. Average execution time of Cetus Parallelized _transf_() routine (with technique applied)
        3. The standard deviation for each of the above mentioned timing results
        4. The performance improvement (Time taken by optimized code/Time taken by the baseline)
        5. Displays the number of threads used


## Translating an input code through the Cetus executable (Optional):

A Cetus executable has been provided which can be used to perform sanity checks ensuring that
Cetus with subscripted subscript analysis enabled generates the expected output code. Note that
only Serial Baseline source codes can be translated.

To perform the sanity check follow the steps below:

```
1. A cetus executable has been provided in the bin subdirectory of
  'Cetus-Subscripted-Subscript_impl' parent directory.
2. Navigate to the source code of the Serial Baseline to be translated:
3. Run the script run-cetus.py:
    $python3 run-cetus.py
4. A new directory 'cetus_output' is created within the same folder as run-cetus.py
    and the translated files are placed inside this new directory.
```
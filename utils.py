
import math,shutil,subprocess
import matplotlib.pylab as plt
import numpy as np

#Plotting routine for Experiment 1
def plot_data_Exp1(benchmark_name,plot_data_dict, plot_title,xlabel,ylabel,save_direc,plt_color):
  
    plot_data_list = plot_data_dict.items()
   
    #plot_data_list = sorted(plot_data_list) 
    x, y = zip(*plot_data_list) 

    plt.plot(x, y,color=plt_color)
    plt.title(plot_title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel) 

    plt.savefig(save_direc + "/" + benchmark_name + '.png')
    plt.close()

#Plotting routine for Experiment 2
def plot_data_Exp2(plot_data_dict:dict, plot_title,save_direc,xlabel,ylabel):

    benchmarks_list = list(plot_data_dict.keys())
    BaseTech_Speedups = []
    NewTech_Speedups = []

    for benchmark in benchmarks_list:
        speedup_tuple = plot_data_dict[benchmark]
        BaseTech_Speedups.append(speedup_tuple[0])
        NewTech_Speedups.append(speedup_tuple[1])
  
    X_axis = np.arange(len(benchmarks_list)) 
  
    plt.figure().set_figwidth(17)
    plt.bar(X_axis - 0.2, BaseTech_Speedups, 0.4, label = 'Cetus+BaseAlgo') 
    plt.bar(X_axis + 0.2, NewTech_Speedups, 0.4, label = 'Cetus+NewAlgo') 
  
    plt.xticks(X_axis, benchmarks_list) 
    plt.ylabel(ylabel) 
    plt.title(plot_title) 
    plt.legend() 

    plt.savefig(save_direc + "Exp2")

#Calculates run-to-run variation (%)
def calculate_variation(values_list, mean):

    sum_of_squares = 0.0

    for val in values_list:
        sum_of_squares += (val - mean)**2 

    standard_dev = math.sqrt(sum_of_squares/len(values_list))
    percent_dev = (100*standard_dev)/mean 
    return percent_dev


#Moves files from a source directory to destination
def move_file(source,destination):
    shutil.move(source,destination)


#Command used to download input matrices
def runcmd(cmd, verbose = False):

    process = subprocess.Popen(
        cmd,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE,
        text = True,
        shell = True
    )
    std_out, std_err = process.communicate()
    # if verbose:
    #     print(std_out.strip(), std_err)
    pass


#Replacing hard coded paths to BLAS and LAPACK binaries
def set_path_to_BLAS_LAPACK(root_path:str,suitesparse_path:str):

    BLAS_PATH= root_path + '/BLAS_and_LAPACK/libopenblas.a -lgfortran'
    LAPACK_PATH = root_path + '/BLAS_and_LAPACK/liblapack.a'

    LAPACK_Replace_String = '~/Artifact-SubSubAnalysis/BLAS_and_LAPACK/liblapack.a'
    BLAS_Replace_String = '~/Artifact-SubSubAnalysis/BLAS_and_LAPACK/libopenblas.a -lgfortran'

    MK_SuiteSparse = suitesparse_path + '/SuiteSparse_config/SuiteSparse_config.mk'

    with open(MK_SuiteSparse, 'r') as file:

        mk_content = file.read()

        if(BLAS_PATH not in mk_content):
            mk_content = mk_content.replace(BLAS_Replace_String, BLAS_PATH) 

        if(LAPACK_PATH not in mk_content):
            mk_content = mk_content.replace(LAPACK_Replace_String, LAPACK_PATH)

    
    with open(MK_SuiteSparse, 'w') as file:

        file.write(mk_content)

        file.close()

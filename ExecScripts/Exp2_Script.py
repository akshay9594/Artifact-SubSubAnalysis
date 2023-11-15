

import matplotlib.pylab as plt

from subprocess import Popen,PIPE

import re,math,os,sys


#Compile the SDDMM benchmark
def compile_fdtd(base_path):
    os.chdir(base_path)
    make_result = Popen('make',stdout=PIPE,stderr=PIPE)
    output, err = make_result.communicate()
    output = str(output, 'UTF-8')
    error = str(err, 'UTF-8')

    if(('error:' in output) or ('error:' in error)):
        print("Compilation failed...")
        exit()
    return



#Running the experiment with amgmk benchmark
def run_exp_fdtd(Exp2_directory,iters,path_to_reports_dir):

    #If the user wants to run the experiment for all the input files

    head_String = "\n===============Timing Results for the fdtd benchmark(Average of "+ str(iters)+" runs)===============\n\n"

    #Set the path of the baseline codes
    base_path = Exp2_directory +'fdtd-2d/Serial/'

    #Set the path of the optimized codes
    opt_code_path = Exp2_directory+'fdtd-2d/Technique_Applied/'


    #Compile the baseline code
    compile_fdtd(base_path)

    #Compile the optimized code
    compile_fdtd(opt_code_path)

    os.chdir(base_path)

    #List of input matrices
    #input_matrices = ['gsm_106857.mtx', 'dielFilterV2clx.mtx','af_shell1.mtx','inline_1.mtx']

   
    with open(path_to_reports_dir+'/SDDMM.txt', 'w') as f:
        f.write(head_String)
        f.write("(a) Baseline Code : Serial code\n")
        f.write("(b) Optimized Code : Cetus Parallel code (with technique applied)\n")

        f.write("\n--------------------------------------------------------------------------\n")

        # For each matrix,execute the baseline code and optimized code and calculate the speedup
    
        #Clean the object files
        os.chdir(base_path)
        Popen(['make','clean'],stdout=PIPE,stderr=PIPE)

        f.write("\n-------------------------------------------------------------------------------\n")
    
    return




#Call specific subroutines that execute the benchmarks
def run_benchmark(benchmark,Exp2_directory,iters,path_to_reports_dir):

    if(benchmark == 'fdtd-2d'):
        run_exp_fdtd(Exp2_directory,iters,path_to_reports_dir)
    else:
        print("Benchmark not supported")
    
    
    return

#Run the main experiment
def RunExp(root_directory):

    iters =  1

    Exp2_directory = root_directory + '/Experiment_2/'

    list_benchmarks = {'poly':['fdtd-2d','heat-3d', 'gramschmidt', 'syrk'],
                        'NAS':['CG', 'MG', 'UA']}
    # list_benchmarks = ['SDDMM']

    for i in range(0,len(list_benchmarks)):

        benchmark = list_benchmarks[i]

        print(str(i+1) + "." , "For Benchmark:", benchmark)

        path_to_reports_dir = os.getcwd() + '/Reports/Experiment_2/' + benchmark

        if(os.path.exists(path_to_reports_dir) == False):
            os.mkdir(path_to_reports_dir)

        run_benchmark(benchmark,Exp2_directory,iters,path_to_reports_dir)

        os.chdir(root_directory)

    print("Experiment finished and Results written to the Reports directory!!")



    return
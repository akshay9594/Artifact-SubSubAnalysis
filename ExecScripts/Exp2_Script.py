

import matplotlib.pylab as plt

from subprocess import Popen,PIPE

import re,math,os,sys


#Compile the SDDMM benchmark
def compile_poly(base_path):
    os.chdir(base_path)
    make_result = Popen('make',stdout=PIPE,stderr=PIPE)
    output, err = make_result.communicate()
    output = str(output, 'UTF-8')
    error = str(err, 'UTF-8')

    if(('error:' in output) or ('error:' in error)):
        print("Compilation failed...")
        exit()
    return


def execute_poly(executable_path,executable):

    exec_result = Popen([executable_path,executable],stdout=PIPE,stderr=PIPE)
    output, err_val = exec_result.communicate()
                
    output = str(output,'UTF-8')
    error = str(err_val, 'UTF-8')
   
    search = re.search('Normalized time: (.*)', output)
    wall_time = [float(i) for i in re.findall("\d+\.\d+", search.group(1))]
    app_time = wall_time.pop()

    search = re.search('Maximal deviation from arithmetic mean of 3 average runs: (.*)', output)
    var = [float(i) for i in re.findall("\d+\.\d+", search.group(1))]
    percent_app_time_var = var.pop()


    return app_time,percent_app_time_var


#Running the experiment with amgmk benchmark
def run_poly_benchmark(Exp2_directory,benchmark,iters,path_to_reports_dir):

    #If the user wants to run the experiment for all the input files

    head_String = "\n===============Timing Results for the fdtd benchmark(Average of "+ str(iters)+" runs)===============\n\n"

    #Set the path of the baseline codes
    base_path = Exp2_directory + benchmark + '/Serial/'

    #Set the path of the optimized codes
    opt_code_path = Exp2_directory + benchmark + '/Techniques_Applied/'

    #Path to benchmark run script
    exec_script = Exp2_directory + 'utilities/./time_benchmark.sh'

    #Path to Serial executable
    serial_executable = base_path + benchmark

    #Path to Opt executable
    opt_executable = opt_code_path + benchmark

    #Compile the baseline code
    compile_poly(base_path)

    #Compile the optimized code
    compile_poly(opt_code_path)

    #Change to base Experiment 2 directory
    os.chdir(Exp2_directory)

    with open(path_to_reports_dir+ '/' + benchmark + '.txt', 'w') as f:
        f.write(head_String)
        f.write("(a) Baseline Code : Serial code\n")
        f.write("(b) Optimized Code : Cetus Parallel code (with technique applied)\n")

        f.write("\n--------------------------------------------------------------------------\n")

        # For each matrix,execute the baseline code and optimized code and calculate the speedup

        app_time,app_time_var = execute_poly(exec_script,serial_executable)

        f.write("->Baseline execution time="+ str(app_time)+" s " + "(" + str(app_time_var)+" % variation)\n")

        opt_app_time,opt_app_time_var = execute_poly(exec_script,opt_executable)

        f.write("->Optimized Code execution time ="+ str(opt_app_time)+" s " + "(" + str(opt_app_time_var)+" % variation)\n")

        app_speedup = app_time/opt_app_time
        f.write("->Speedup="+str(app_speedup)+"\n")
        f.write("\n-------------------------------------------------------------------------------\n")

    #Clean the object files
    os.chdir(base_path)
    Popen(['make','clean'],stdout=PIPE,stderr=PIPE)

    os.chdir(opt_code_path)
    Popen(['make','clean'],stdout=PIPE,stderr=PIPE)
    
    return




#Call specific subroutines that execute the benchmarks
def drive_poly(Exp2_directory,root_directory,iters,list_benchmarks):

    for i in range(0,len(list_benchmarks)):
        
        benchmark = list_benchmarks[i]

        print(str(i+1) + "." , "For Benchmark:", benchmark)

        path_to_reports_dir = os.getcwd() + '/Reports/Experiment_2/' + benchmark

        if(os.path.exists(path_to_reports_dir) == False):
            os.mkdir(path_to_reports_dir)

        run_poly_benchmark(Exp2_directory,benchmark,iters,path_to_reports_dir)

        os.chdir(root_directory)

    return

#Run the main experiment
def RunExp(root_directory):

    iters =  1

    Exp2_directory = root_directory + '/Experiment_2/'

    benchmarks_dict = {'poly':['fdtd-2d','heat-3d', 'gramschmidt', 'syrk'],
                        'NAS':['CG', 'MG', 'UA', 'IS'], 'rest':['SDDMM','ic0_csc']}
    
    # list_benchmarks = ['SDDMM']

    benchmark_tags = list(benchmarks_dict.keys())

    for tag in benchmark_tags:
        if(tag == 'poly'):
            list_benchmarks = benchmarks_dict[tag]
            drive_poly(Exp2_directory,root_directory,iters,list_benchmarks)


    print("Experiment finished and Results written to the Reports directory!!")



    return
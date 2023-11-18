

import matplotlib.pylab as plt

from subprocess import Popen,PIPE

import re,math,os,sys

import utils

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

def compile_NAS(code_type:str,base_path:str,benchmark:str,input_class:str):
    os.chdir(base_path)

    make_result = ''
   
    if(code_type == 'baseline'):
        make_result = Popen(['make',benchmark,'CLASS='+input_class],stdout=PIPE,stderr=PIPE)
    elif(code_type == 'opt'):
        make_result = Popen('make',stdout=PIPE,stderr=PIPE)
    else:
        make_result = Popen('make',stdout=PIPE,stderr=PIPE)
    
    output, err = make_result.communicate()
    output = str(output, 'UTF-8')
    error = str(err, 'UTF-8')

    if(('error:' in output) or ('error:' in error)):
        print("Compilation failed...")
        print(output)
        print(error)
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


#Executing the UA application
def execute_NAS(exec_path,exec_command,benchmark,iters):
    transf_times = []
    os.chdir(exec_path)

    search_string = ''
    if(benchmark == 'UA'):
        search_string = 'Transf total Time=(.*)seconds'
    else:
        search_string = 'Time in seconds =(.*)'

    for i in range(0,iters):
        exec_result = Popen(exec_command,stdout=PIPE,stderr=PIPE)
        output, err_val = exec_result.communicate()
                
        output = str(output,'UTF-8')
       
        if(isinstance(output,str)):
            search = re.search(search_string, output)
            wall_time = [float(i) for i in re.findall("\d+\.\d+", search.group(1))]
            transf_times.append(wall_time.pop())

            
    transf_time_sum = 0.0
    for i in range(0,len(transf_times)):
        transf_time_sum += transf_times[i]
   
    mean = transf_time_sum/iters

    percent_var = utils.calculate_variation(transf_times,mean)

    return (mean,percent_var)


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


#Running the experiment with amgmk benchmark
def run_NAS_benchmark(Exp2_directory,benchmark,iters,path_to_reports_dir,cl,executable):

    #If the user wants to run the experiment for all the input files

    head_String = "\n===============Timing Results for the Kernel transf NAS-UA benchmark(Average of "+ str(iters)+" runs)===============\n\n"

    with open(path_to_reports_dir+'/NAS-'+ benchmark +'.txt', 'w') as f:
        f.write(head_String)
        f.write("(a) Baseline Code : Serial kernel\n")
        f.write("(b) Optimized Code : Cetus Parallel kernel (with technique applied)\n")

        f.write("\n--------------------------------------------------------------------------\n")

        
        f.write("Input Class: "+cl+"\n")

        #Set the path of the baseline codes
        base_path = Exp2_directory +'/'+benchmark+'-NAS/Serial/'

        #Compile the baseline code
        compile_NAS('baseline',base_path,benchmark,cl)

        #Clean the object files
        os.chdir(base_path)
        Popen(['make','clean'],stdout=PIPE,stderr=PIPE)

        #Execute the baseline code
        exec_path = base_path + 'bin/'
        exec_command = './'+executable
        app_time, app_time_var = execute_NAS(exec_path,exec_command,benchmark,iters)

        f.write("->Baseline execution time ="+ str(app_time)+" s " + "(" + str(app_time_var)+" % variation)\n")
           

        #Set the path of the optimized codes
        opt_code_path = Exp2_directory+'/'+ benchmark +'-NAS/Techniques_Applied/'

        #Compile and execute the optimized code
        compile_NAS('opt',opt_code_path,benchmark,cl)
        opt_app_time, opt_time_var = execute_NAS(opt_code_path,exec_command,benchmark,iters)

        f.write("->Optimized Code execution time ="+ str(opt_app_time)+" s " + "(" + str(opt_time_var)+" % variation)\n")

        app_speedup = app_time/opt_app_time
        f.write("->Speedup="+str(app_speedup)+"\n")

        #Clean the object files
        os.chdir(base_path)
        Popen(['make','clean'],stdout=PIPE,stderr=PIPE)

        os.chdir(opt_code_path)
        Popen(['make','clean'],stdout=PIPE,stderr=PIPE)


        f.write("\n-------------------------------------------------------------------------------\n")

    return


#Driver code to execute the Polybench benchmarks
def drive_poly(Exp2_directory,root_directory,iters,list_benchmarks):

    for i in range(0,len(list_benchmarks)):
        
        benchmark = list_benchmarks[i]

        print(str(i+1) + "." , "For Benchmark:", benchmark)

        path_to_reports_dir = os.getcwd() + '/Reports/Experiment_2/'

        if(os.path.exists(path_to_reports_dir) == False):
            os.mkdir(path_to_reports_dir)

        #Actual subroutine that executes the benchmark
        run_poly_benchmark(Exp2_directory,benchmark,iters,path_to_reports_dir)

        os.chdir(root_directory)

    return

#Driver code that executes the NAS benchmarks
def drive_NAS(Exp2_directory,root_directory,iters,list_benchmarks):

    for i in range(0,len(list_benchmarks)):
        
        benchmark = list_benchmarks[i]

        print(str(i+1) + "." , "For Benchmark:", benchmark)

        path_to_reports_dir = os.getcwd() + '/Reports/Experiment_2/'

        if(os.path.exists(path_to_reports_dir) == False):
            os.mkdir(path_to_reports_dir)

    #Set the input class to be used and executable for each benchmark
        input_class = ''
        baseline_executable = ''
        if(benchmark == 'UA'):
            input_class = 'A'
            executable = 'ua.' + input_class + '.x'
        elif(benchmark == 'MG'):
            input_class = 'B'
            executable = 'mg.' + input_class + '.x'
        elif(benchmark == 'CG'):
            input_class = 'B'
            executable = 'cg.' + input_class + '.x'
        elif(benchmark == 'IS'):
            input_class = 'C'
            executable = 'is.' + input_class + '.x'
        else:
            print("NAS benchmark not supported!!")
            return

        #Actual subroutine that executes the benchmark
        run_NAS_benchmark(Exp2_directory,benchmark,iters,path_to_reports_dir,input_class,executable)

        os.chdir(root_directory)

    return


#Run the main experiment
def RunExp(root_directory):

    iters =  1

    Exp2_directory = root_directory + '/Experiment_2/'

    #Benchmarks in a dictionary. They are grouped according to the suite they belong.
    #Grouping helps in seamleass compilation and execution.
    benchmarks_dict = {'poly':['fdtd-2d','heat-3d', 'gramschmidt', 'syrk'],
                        'NAS':['CG', 'MG', 'UA', 'IS'], 'rest':['SDDMM','ic0_csc']}
    

    benchmark_tags = list(benchmarks_dict.keys())

    for tag in benchmark_tags:
        list_benchmarks = benchmarks_dict[tag]
        # if(tag == 'poly'):
        #     #drive_poly(Exp2_directory,root_directory,iters,list_benchmarks)

        if(tag == 'NAS'):
            drive_NAS(Exp2_directory,root_directory,iters,list_benchmarks)


    print("Experiment finished and Results written to the Reports directory!!")



    return
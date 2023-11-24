
import matplotlib.pylab as plt

from subprocess import Popen,PIPE

import re,os,subprocess

import utils

#Executing the AMGmk application
def execute_amgmk(iters):
    app_times = []

    threads = 0
    for i in range(0,iters):
        exec_result = Popen('./AMGMk',stdout=PIPE,stderr=PIPE)
        output, err_val = exec_result.communicate()
                
        output = str(output,'UTF-8')
        error = str(err_val, 'UTF-8')
        if(isinstance(output,str)):
            search = re.search('Wall time =(.*)seconds.', output)
            wall_time = [float(i) for i in re.findall("\d+\.\d+", search.group(1))]
            app_times.append(wall_time.pop())


            if(i == 0):
                search = re.search('max_num_threads =(.*)', output)
                threads = [int(s) for s in re.findall(r'\b\d+\b', search.group(1))].pop()
            
    app_time_sum = 0.0
    for i in range(0,len(app_times)):
        app_time_sum += app_times[i]

    app_time_avg = app_time_sum/iters

    percent_app_time_var = utils.calculate_variation(app_times,app_time_avg)

    return (app_time_avg,percent_app_time_var,threads)


#Executing the UA application
def execute_UA(exec_path,input_class,iters):
    transf_times = []
    os.chdir(exec_path)

    exec_command = './' + 'ua.'+input_class+'.x'
    threads = 0
    for i in range(0,iters):
        exec_result = Popen(exec_command,stdout=PIPE,stderr=PIPE)
        output, err_val = exec_result.communicate()
                
        output = str(output,'UTF-8')
       
        if(isinstance(output,str)):
            search = re.search('Transf total Time=(.*)seconds', output)
            wall_time = [float(i) for i in re.findall("\d+\.\d+", search.group(1))]
            transf_times.append(wall_time.pop())

            if(i == 0):
                search = re.search('max_num_threads =(.*)', output)
                threads = [int(s) for s in re.findall(r'\b\d+\b', search.group(1))].pop()
            
    transf_time_sum = 0.0
    for i in range(0,len(transf_times)):
        transf_time_sum += transf_times[i]
   
    mean = transf_time_sum/iters

    percent_var = utils.calculate_variation(transf_times,mean)

    return (mean,percent_var,threads)


#Executing the SDDMM application
def execute_SDDMM(executable_path,input_path,iters):
    app_times = []

    threads = 0

    executable_path = executable_path + './sddmm'
   
    for i in range(0,iters):
        exec_result = Popen([executable_path,input_path],stdout=PIPE,stderr=PIPE)
        output, err_val = exec_result.communicate()
                
        output = str(output,'UTF-8')
        error = str(err_val, 'UTF-8')
        if(isinstance(output,str)):
            search = re.search('kernel=(.*) s', output)
            wall_time = [float(i) for i in re.findall("\d+\.\d+", search.group(1))]
            app_times.append(wall_time.pop())


            if(i == 0):
                search = re.search('max_num_threads =(.*)', output)
                threads = [int(s) for s in re.findall(r'\b\d+\b', search.group(1))].pop()
            
    app_time_sum = 0.0
    for i in range(0,len(app_times)):
        app_time_sum += app_times[i]

    app_time_avg = app_time_sum/iters

    percent_app_time_var = utils.calculate_variation(app_times,app_time_avg)

    return (app_time_avg,percent_app_time_var,threads)



#Compile the AMGmk benchmark
def compile_amgmk(base_path):
    os.chdir(base_path)
    make_result = Popen('make',stdout=PIPE,stderr=PIPE)
    output, err = make_result.communicate()
    output = str(output, 'UTF-8')
    error = str(err, 'UTF-8')

    if(('error:' in output) or ('error:' in error)):
        print("Compilation failed...")
        exit()
    return

#Compile the UA benchmark
def compile_UA(code_type,base_path,input_class):
    os.chdir(base_path)

    make_result = ''
   
    if(code_type == 'baseline'):
        make_result = Popen(['make','UA','CLASS='+input_class],stdout=PIPE,stderr=PIPE)
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

#Compile the SDDMM benchmark
def compile_SDDMM(base_path):
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
def run_exp_amgmk(Exp1_directory,iters,path_to_reports_dir):

    #If the user wants to run the experiment for all the input files

    head_String = "\n===============Timing Results for the AMGmk benchmark(Average of "+ str(iters)+" runs)===============\n\n"

    with open(path_to_reports_dir+'/AMGmk.txt', 'w') as f:
        f.write(head_String)
        f.write("(a) Baseline Code : Serial Code\n")
        f.write("(b) Optimized Code : Cetus Parallel Code (with technique applied)\n")

        f.write("\n--------------------------------------------------------------------------\n")

        speedup_dict = {}
        for i in range(1,6):
            input_matrix = 'MATRIX'+ str(i)
            #Set the path of the baseline codes
            base_path = Exp1_directory +'/amgmk-v1.0/Baselines/Serial/' + input_matrix

            #Set the path of the optimized codes
            opt_code_path = Exp1_directory+'/amgmk-v1.0/Technique_Applied/' + input_matrix

            f.write(str(i) + ". Input Matrix: " + input_matrix + "\n")

            #Compile the baseline code
            compile_amgmk(base_path)
            #Execute the baseline code
            app_time,app_time_var,threads = execute_amgmk(iters)
            
            #Compile the optimized code
            compile_amgmk(opt_code_path)
            # #Execute the optimized code
            opt_app_time,opt_app_time_var,threads = execute_amgmk(iters)

            f.write("->Baseline execution time="+ str(app_time)+" s " + "(" + str(app_time_var)+" % variation)\n")
           
            f.write("->Optimized Code execution time ="+ str(opt_app_time)+" s " + "(" + str(opt_app_time_var)+" % variation)\n")

            app_speedup = app_time/opt_app_time
            f.write("->Speedup="+str(app_speedup)+"\n")

            speedup_dict[input_matrix] = app_speedup

            #Clean the object files
            os.chdir(base_path)
            Popen(['make','clean'],stdout=PIPE,stderr=PIPE)

            os.chdir(opt_code_path)
            Popen(['make','clean'],stdout=PIPE,stderr=PIPE)


            f.write("\n-------------------------------------------------------------------------------\n")

    return speedup_dict



#Running the experiment with the UA benchmark
def run_exp_UA(Exp1_directory,iters,path_to_reports_dir):

    #If the user wants to run the experiment for all the input files

    head_String = "\n===============Timing Results for the Kernel transf NAS-UA benchmark(Average of "+ str(iters)+" runs)===============\n\n"

    input_classes = ['A', 'B']      # Include Class C

    speedup_dict = {}

    with open(path_to_reports_dir+'/NAS-UA.txt', 'w') as f:
        f.write(head_String)
        f.write("(a) Baseline Code : Serial transf kernel\n")
        f.write("(b) Optimized Code : Cetus Parallel transf kernel (with technique applied)\n")

        f.write("\n--------------------------------------------------------------------------\n")


        for i in range(0,len(input_classes)):
            cl = input_classes[i]
            f.write(str(i+1)+". Input Class: "+cl+"\n")

            #Set the path of the baseline codes
            base_path = Exp1_directory +'/UA-NAS/Baselines/Serial/'

            #Compile the baseline code
            compile_UA('baseline',base_path,cl)

            #Clean the object files
            os.chdir(base_path)
            Popen(['make','clean'],stdout=PIPE,stderr=PIPE)

            #Execute the baseline code
            exec_path = base_path + 'bin/'
            app_time, app_time_var, threads = execute_UA(exec_path,cl,iters)

            f.write("->Baseline execution time ="+ str(app_time)+" s " + "(" + str(app_time_var)+" % variation)\n")
           

            #Set the path of the optimized codes
            opt_code_path = Exp1_directory+'/UA-NAS/Technique_Applied/CLASS-' + cl

            #Compile and execute the optimized code
            compile_UA('opt',opt_code_path,cl)
            opt_app_time, opt_time_var, threads = execute_UA(opt_code_path,cl,iters)

            f.write("->Optimized Code execution time ="+ str(opt_app_time)+" s " + "(" + str(opt_time_var)+" % variation)\n")

            app_speedup = app_time/opt_app_time
            f.write("->Speedup="+str(app_speedup)+"\n")

            speedup_dict['CLASS-'+cl] = app_speedup
            #Clean the object files
            os.chdir(base_path)
            Popen(['make','clean'],stdout=PIPE,stderr=PIPE)

            os.chdir(opt_code_path)
            Popen(['make','clean'],stdout=PIPE,stderr=PIPE)


            f.write("\n-------------------------------------------------------------------------------\n")

    return speedup_dict



#Running the experiment with the SDDMM benchmark
def run_exp_SDDMM(Exp1_directory,iters,path_to_reports_dir):

    #If the user wants to run the experiment for all the input files

    head_String = "\n===============Timing Results for the SDDMM benchmark(Average of "+ str(iters)+" runs)===============\n\n"

    #Set the path of the baseline codes
    base_path = Exp1_directory +'SDDMM/Baselines/Serial/'

    #Set the path of the optimized codes
    opt_code_path = Exp1_directory+'SDDMM/Technique_Applied/'

    input_directory = os.getcwd() +'/input_matrices/'


    #Compile the baseline code
    compile_SDDMM(base_path)

    #Compile the optimized code
    compile_SDDMM(opt_code_path)

    os.chdir(base_path)

    #List of input matrices
    input_matrices = ['gsm_106857', 'dielFilterV2clx','af_shell1','inline_1']

    # input_matrices = ['af_shell1'] 
    speedup_dict = {} 

    with open(path_to_reports_dir+'/SDDMM.txt', 'w') as f:
        f.write(head_String)
        f.write("(a) Baseline Code : Serial code\n")
        f.write("(b) Optimized Code : Cetus Parallel code (with technique applied)\n")

        f.write("\n--------------------------------------------------------------------------\n")

        # For each matrix,execute the baseline code and optimized code and calculate the speedup
        for i in range(0,len(input_matrices)):

            matrix = input_matrices[i]

            f.write(str(i+1)+". For matrix: "+ matrix+"\n")

            #Path to the input matrix
            input_path = input_directory + matrix + "/" + matrix + ".mtx"
            
            #Execute the Baseline code
            app_time, app_time_var, threads = execute_SDDMM(base_path,input_path,iters)
            #Write the execution times to the output file
            f.write("->Baseline execution time ="+ str(app_time)+" s " + "(" + str(app_time_var)+" % variation)" +"\n")
            #Execute the optimized code
            opt_app_time, opt_time_var, threads = execute_SDDMM(opt_code_path,input_path,iters)
            #Writ the execution time to the output file
            f.write("->Optimized Code execution time ="+ str(opt_app_time)+" s " + "(" + str(opt_time_var)+" % variation)"+"\n")
            #Calculate the speedup
            app_speedup = app_time/opt_app_time
            f.write("->Speedup="+str(app_speedup)+"\n")

            speedup_dict[matrix] = app_speedup

            f.write("\n-------------------------------------------------------------------------------\n")

        #Clean the object files
        os.chdir(base_path)
        Popen(['make','clean'],stdout=PIPE,stderr=PIPE)

        os.chdir(opt_code_path)
        Popen(['make','clean'],stdout=PIPE,stderr=PIPE)

    
    return speedup_dict



#Call specific subroutines that execute the benchmarks
def run_benchmark(benchmark,Exp1_directory,iters,path_to_reports_dir):

    speedup_dict = {}
    xlabel = ''
    plot_title = ' '
    if(benchmark == 'AMGmk'):
        speedup_dict = run_exp_amgmk(Exp1_directory,iters,path_to_reports_dir)
        plot_title = "Performance improvement of the AMGMk benchmark"
        xlabel = "Input Matrices"

    elif(benchmark == 'UA-NAS'):
        
        speedup_dict = run_exp_UA(Exp1_directory,iters,path_to_reports_dir)
        plot_title = "Performance improvement of the UA benchmark"
        xlabel = "Input Class"

    elif(benchmark == 'SDDMM'):
        speedup_dict = run_exp_SDDMM(Exp1_directory,iters,path_to_reports_dir)
        xlabel = "Input Matrices"
        plot_title = "Performance improvement of the SDDMM benchmark"
        

    else:
        print("Benchmark not supported")
    
    
    return speedup_dict,plot_title,xlabel


#Run the main experiment
def RunExp(root_directory):

    iters =  1

    Exp1_directory = root_directory + '/Experiment_1/'

    list_benchmarks = ['AMGmk','UA-NAS', 'SDDMM']
    #list_benchmarks = ['SDDMM']

    ylabel="Performance Improvement"
    for i in range(0,len(list_benchmarks)):

        benchmark = list_benchmarks[i]

        print(str(i+1) + "." , "For Benchmark:", benchmark)

        path_to_reports_dir = os.getcwd() + '/Reports/Experiment_1/' + benchmark

        if(os.path.exists(path_to_reports_dir) == False):
            os.mkdir(path_to_reports_dir)

        speedup_dict,plot_title,xlabel = run_benchmark(benchmark,Exp1_directory,iters,path_to_reports_dir)

        os.chdir(root_directory)

        #Plot the speedup data
        utils.plot_data(benchmark,speedup_dict,plot_title,xlabel,ylabel,path_to_reports_dir)

    print("Experiment finished and Results written to the Reports directory!!")



    return
import matplotlib.pylab as plt

from subprocess import Popen,PIPE

import re,math

import os


print("\n*************************************************************************\n")
print("\t\t>Artifact execution for reproducing the results<")
print("\n*************************************************************************\n")



def plot_data(plot_data_dict, plot_title,xlabel,ylabel):
    plot_data_list = plot_data_dict.items()
    plot_data_list = sorted(plot_data_list) 
    x, y = zip(*plot_data_list) 

    plt.bar(x, y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(plot_title)          
    plt.show()


def calculate_variation(values_list, mean):

    sum_of_squares = 0.0

    for val in values_list:
        sum_of_squares += (val - mean)**2 

    standard_dev = math.sqrt(sum_of_squares/len(values_list))
    percent_dev = (100*standard_dev)/mean 
    return percent_dev

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


def compile_UA(baseline,base_path,input_class):
    os.chdir(base_path)

    make_result = ''
   
    if(baseline == '1'):
        make_result = Popen(['make','CLASS='+input_class],stdout=PIPE,stderr=PIPE)
    elif(baseline == '2'):
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

#Executing the amgmk application
def execute_amgmk(iters):
    app_times = []
    loop_times = []

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

            search = re.search('Target loop time=(.*)seconds', output)
            loop_time = [float(i) for i in re.findall("\d+\.\d+", search.group(1))]
            loop_times.append(loop_time.pop())

            if(i == 0):
                search = re.search('max_num_threads =(.*)', output)
                threads = [int(s) for s in re.findall(r'\b\d+\b', search.group(1))].pop()
            
    app_time_sum = 0.0
    for i in range(0,len(app_times)):
        app_time_sum += app_times[i]
    loop_time_sum = 0.0
    for i in range(0,len(loop_times)):
        loop_time_sum += loop_times[i]

    app_time_avg = app_time_sum/iters
    loop_time_avg = loop_time_sum/iters

    percent_app_time_var = calculate_variation(app_times,app_time_avg)
    percent_loop_time_var = calculate_variation(loop_times,loop_time_avg)

    Popen(['make','clean'],stdout=PIPE,stderr=PIPE)
    return (app_time_avg,percent_app_time_var,loop_time_avg,percent_loop_time_var,threads)

#Executing the UA application
def execute_UA(exec_path,input_class,iters):
    transf_times = []
    os.chdir(exec_path)

    exec_command = './ua.'+input_class+'.x'
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

    percent_var = calculate_variation(transf_times,mean)

    return (mean,percent_var,threads)


#Selecting the base line:
#There are two possible baselines -
#(1) The serial code
#(2) The Cetus output code (Without subscripted subscript analysis of the technique applied)
def select_baseline():

    baseline = input("\n\t1.Serial\n\t2.Cetus-Output-WithoutSubSub (Takes substantial run-time)\n\t->Select the baseline for comparison:\n")
    if(baseline == '1'):
        print("{Serial baseline selected}\n")
    elif(baseline == '2'):
        print("{Cetus Output without technique applied selected as baseline}\n")
    else:
        print('Invalid Input')
        exit()

    return baseline


#Running the experiment with amgmk benchmark
def run_exp_amgmk(base_path,opt_code_path,iters,f):

    #Compile the baseline code
    compile_amgmk(base_path)
    #Execute the baseline code
    app_time,app_time_var,loop_time,loop_time_var,threads = execute_amgmk(iters)
      
    #Compile the optimized code
    compile_amgmk(opt_code_path)
    #Execute the optimized code
    opt_app_time,opt_app_time_var,opt_loop_time,opt_loop_time_var,threads = execute_amgmk(iters)

    #Write the absolute numbers to the output text file
    f.write("(a) For the Application:\n")
    f.write("->Average baseline application execution time="+ str(app_time)+" s\n")
    f.write("->Baseline application execution time Variation="+ str(app_time_var)+" %\n\n")
    f.write("->Average application execution time of Cetus Parallelized Code (with technique applied)="+ str(opt_app_time)+ " s\n")
    f.write("->Optimized application code execution time Variation="+ str(opt_app_time_var)+" %\n")
    app_speedup = app_time/opt_app_time
    f.write("->Application Speedup="+str(app_speedup)+"\n")
    
    f.write("\n(b) For the kernel:\n")
    f.write("->Average baseline kernel execution time="+ str(loop_time)+ " s\n")
    f.write("->Baseline kernel execution time Variation="+ str(loop_time_var)+" %\n\n")    
    f.write("->Average kernel execution time of Cetus Parallelized Code (with technique applied)="+str(opt_loop_time)+" s\n")
    f.write("->Optimized kernel code execution time Variation="+ str(opt_loop_time_var)+" %\n")
    knl_speedup = loop_time/opt_loop_time
    f.write("->Kernel Speedup=" + str(knl_speedup)+"\n")
    f.write("->threads used=" + str(threads)+"\n")

    return app_speedup,knl_speedup,threads


#Running the experiment with the UA benchmark
def run_exp_UA(baseline,base_path,opt_code_path,input_class,iters,f):

    #Change pth to the appropriate directory which has the baseline code
    if(baseline == '1'):
        compile_path = base_path + 'UA'
        exec_path = base_path +'bin'
    elif(baseline == '2'):
        compile_path = base_path + 'CLASS-'+input_class
        exec_path = compile_path

    # #Compile and execute the baseline code
    compile_UA(baseline,compile_path,input_class)
    avg_base_time,variation,threads = execute_UA(exec_path,input_class,iters)

    f.write("->Average baseline subroutine execution time="+str(avg_base_time)+ " s\n")
    f.write("->Baseline Execution time variation="+str(variation)+ " %\n")
    
    #Clean the baseline object files
    os.chdir(base_path)
    Popen(['make','clean'],stdout=PIPE,stderr=PIPE)

    #Change to the directory which has the optimized code (Parallel code with technique applied)
    compile_path = opt_code_path
    exec_path = compile_path

    #Compile and execute the code
    baseline = 'NULL'
    compile_UA(baseline,compile_path,input_class)
    avg_opt_time,variation,threads = execute_UA(exec_path,input_class,iters) 
    
    subroutine_speedup = avg_base_time/avg_opt_time

    f.write("->Average subroutine execution time of Cetus Parallelized Code (with technique applied) = "+str(avg_opt_time)+" s\n")
    f.write("->Optimized Code Execution time variation="+str(variation)+ " %\n")
    f.write("->Subroutine Speedup = "+str(subroutine_speedup)+"\n")
    f.write("->threads used = "+str(threads)+"\n")

    #Clean the optimized code object files
    os.chdir(compile_path)
    Popen(['make','clean'],stdout=PIPE,stderr=PIPE)

    return subroutine_speedup,threads






#Main code
#Step 1: Select the benchmark to be evaluated
val = input("\t1. amgmk-v1.0\n\t2. UA-NAS\n\t->Select the benchmark you want to test by entering the number:\n")

#If amgmk-v1.0 selected for evaluation
if(val == '1'):
    print("{amgmk-v1.0 selected for evaluation}\n")
    input_matrix = ''
    #Step 2: Choose the baseline
    baseline = select_baseline()
    base_path = ''
    root = os.getcwd()

    #Set the path to the Baseline source codes
    if(baseline == '1'):
        base_path = root+'/amgmk-v1.0/Baselines/Serial/'
            
    elif(baseline == '2'):
        base_path = root+'/amgmk-v1.0/Baselines/Cetus-Output-WithoutSubSub/'

    opt_code_path = root+'/amgmk-v1.0/Technique_Applied/'

    application_speedups = {}
    kernel_speedups = {}

    #Step 3: Choose if you want to run 1 matrix or all matrices
    ans = input("\t->Do you want to run the experiment for all input matrices?yes or no\n")
    if(ans != 'yes' and ans != 'no'):
        print("invalid input")
        exit()

    #Step 4: Specify the number of application runs
    iters = int(input("\t->Specify the number of application runs:\n"))
    head_String = "\n=========================Timing Results for the amgmk benchmark(Average of "+ str(iters)+" runs)===============================\n\n"
    
    #Open a text file to write the results
    with open('output.txt', 'w') as f:
        f.write(head_String)
        #If the user wants to run the experiment for all the input files
        if(ans == 'yes'):
            for i in range(1,6):
                input_matrix = 'MATRIX'+ str(i)

                #Set the path of the baseline codes
                compile_path = base_path + input_matrix
                f.write(str(i)+". For Input matrix: "+input_matrix+"\n")

                #Set the path of the optimized codes
                Technique_Applied_path = opt_code_path + input_matrix

                #Run the experiment and gather the application and kernel speedups
                app_speedup,knl_speedup,num_threads = run_exp_amgmk(compile_path,Technique_Applied_path,iters,f)
                application_speedups[input_matrix] = app_speedup
                kernel_speedups[input_matrix] = knl_speedup
                f.write("------------------------------------------------------------------------------------------------------\n")
          
           #Plot the Speedup data for the input matrices
            if(baseline == '1'):
                plt_app_title = 'Performance improvement of the parallel application (Technique Applied v/s Serial) on ' +str(num_threads)+' threads'
                plt_kernel_title = 'Performance improvement of the parallel kernel (Technique Applied v/s Serial) on ' +str(num_threads)+' threads'
            else:
                plt_app_title = 'Performance improvement of the parallel application (Technique Applied v/s Technique NOT Applied) on ' +str(num_threads)+' threads'
                plt_kernel_title = 'Performance improvement of the parallel kernel (Technique Applied v/s Technique NOT Applied) on ' +str(num_threads)+' threads'
            
            xlabel = 'Input Matrices'
            ylabel = 'Performance Improvement'

            plot_data(application_speedups,plt_app_title,xlabel,ylabel)
            plot_data(kernel_speedups,plt_kernel_title,xlabel,ylabel)

        #If the user wants to run the experiment for a specific input matrix
        #The speedup in this case is not plotted
        elif(ans == 'no'):
            mat_num = input("\n->Enter the number of the matrix from 1 to 5 to use as input:")
            input_matrix = 'MATRIX' + str(mat_num)
            print("{Input matrix selected for experiment: ",input_matrix,"}\n")

            base_path += input_matrix
            opt_code_path += input_matrix
            app_speedup,kernel_speedup,num_threads = run_exp_amgmk(base_path,opt_code_path,iters,f)

            application_speedups[input_matrix] = app_speedup
            kernel_speedups[input_matrix] = kernel_speedup

        else:
            print("Invalid Input")
            exit()

        os.chdir(root)
        print("Experiment completed successfully! Results have been written to output.txt")

elif(val == '2'):
    print("{UA-NAS selected for evaluation}")

    input_class = ''
    baseline = select_baseline()
    ans = input("\t->Do you want to run the experiment for all input classes?yes or no\n")
    if(ans != 'yes' and ans != 'no'):
        print("invalid input")
        exit()

    base_path = ''
    root = os.getcwd()

    if(baseline == '1'):
        base_path = root+'/UA-NAS/Baselines/Serial/'
            
    elif(baseline == '2'):
        base_path = root+'/UA-NAS/Baselines/Cetus-Output-WithoutSubSub/'
    
    opt_code_path = root+'/UA-NAS/Technique_Applied/'

    iters = int(input("\t->Specify the number of application runs:\n"))
    head_String = "\n=========================Timing Results for the transf routine from the UA benchmark(Average of "+ str(iters)+" runs)===============================\n\n"

    with open('output.txt', 'w') as f:
        f.write(head_String)
        f.write('\t***Note: Input Class D not included since it takes about 2 to 3 hours to complete***\n\n')

        if(ans == 'yes'):
            input_classes = ['A', 'B', 'C']
            count = 1
            speedups = {}
            for cl in input_classes:
                f.write(str(count)+". For Input Class: "+cl+"\n")
                opt_code_for_class = opt_code_path + 'CLASS-'+cl
                subroutine_speedup,num_threads = run_exp_UA(baseline,base_path,opt_code_for_class,cl,iters,f)
                speedups['Class '+cl] = subroutine_speedup
                count = count + 1
                f.write("------------------------------------------------------------------------------------------------------\n")

         #Plot the Speedup data for the input matrices
            if(baseline == '1'):
                plt_title = 'Performance improvement of the transf routine (Technique Applied v/s Serial) on '+str(num_threads)+' threads'
            else:
                plt_title = 'Performance improvement of the transf routine (Technique Applied v/s Technique NOT Applied) on '+str(num_threads)+' threads'

            xlabel = 'Input Classes'
            ylabel = 'Performance Improvement'

            plot_data(speedups,plt_title,xlabel,ylabel)

        elif(ans == 'no'):
            input_class = input("Enter the Input Class (A,B,C or D) to use:")
            opt_code_path += 'CLASS-'+input_class

            subroutine_speedup = run_exp_UA(baseline,base_path,opt_code_path,input_class,iters,f)

        else:
            print("Invalid Input")
            exit()

    os.chdir(root)
    print("Experiment completed successfully! Results have been written to output.txt")

else:
    print("Invalid input")


import math,shutil,subprocess
import matplotlib.pylab as plt
import os

def plot_data(benchmark_name,plot_data_dict, plot_title,xlabel,ylabel,save_direc):
    plot_data_list = plot_data_dict.items()
    #plot_data_list = sorted(plot_data_list) 
    x, y = zip(*plot_data_list) 

    plt.plot(x, y)
    plt.title(plot_title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel) 

    plt.savefig(save_direc + "/" + benchmark_name + '.png')

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
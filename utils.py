
import math

#Calculates run-to-run variation (%)
def calculate_variation(values_list, mean):

    sum_of_squares = 0.0

    for val in values_list:
        sum_of_squares += (val - mean)**2 

    standard_dev = math.sqrt(sum_of_squares/len(values_list))
    percent_dev = (100*standard_dev)/mean 
    return percent_dev
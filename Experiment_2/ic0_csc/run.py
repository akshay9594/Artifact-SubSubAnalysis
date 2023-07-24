

import subprocess, sys


print("Testing the incomplete Cholesky Factorization Benchmark")


val = input("Compile and Run the serial code or the Cetus parallel code?\n1. Serial\n2. Cetus Parallel\n ")

test_files = ''
if(val == '1'):
    test_files = 'IncompleteCholesky_m7.c'
    
elif(val == '2'):
    test_files = './cetus_output/IncompleteCholesky_m7.c'
   
else:
    print("Invalid Input (Type either 1 or 2")
    sys.exit()


test_result = subprocess.call(['gcc', '-fopenmp' ,'-o', 'icm7', test_files, '-lm'])#stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
 

if(test_result != 0):
    print("Compilation failed!!")

else:
    print("Compilation Succeeded...")

input = 'bmwcra_1.mtx'

print("Testing with input file: bmwcra_1.mtx:")

test_result = subprocess.call(['./icm7', 'input/'+input])#stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
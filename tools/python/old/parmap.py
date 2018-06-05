from multiprocessing import Process, Pipe, cpu_count 
from itertools import izip  
import numpy as np

def spawn(f):  
    def fun(pipe,x):  
        pipe.send(f(x))  
        pipe.close()  
    return fun  

def parmap(f,X,np=None):  
    if np == None:
        np = cpu_count()/2

    outputList = []  
    pipe=[Pipe() for x in X]  
    processes=[Process(target=spawn(f),args=(c,x)) for x,(p,c) in izip(X,pipe)]  

    numProcesses = len(processes)  
    processNum = 0  
    while processNum < numProcesses:  
        endProcessNum = min(processNum+np, numProcesses)  
        for proc in processes[processNum:endProcessNum]:  
            proc.start()  
        for proc in processes[processNum:endProcessNum]:  
            proc.join(1.)  
        for proc,c in pipe[processNum:endProcessNum]:  
            outputList.append(proc.recv())  

        processNum = endProcessNum  
                    
    return outputList 


def parmapLimitThreads(f,X,np=None,maxThreads=120):  
    if np == None:
        np = cpu_count()/2

    numThreads = len(X)
    nthreads = min(numThreads,maxThreads)
    threadNum = 0
    outputList = []  
    while threadNum < numThreads:
        endThreadNum = min(threadNum+nthreads, numThreads)
        pipe=[Pipe() for x in X[threadNum:endThreadNum]]  
        processes=[Process(target=spawn(f),args=(c,x)) for x,(p,c) in izip(X[threadNum:endThreadNum],pipe)]  
    
        numProcesses = len(processes)  
        processNum = 0  
        while processNum < numProcesses:  
            endProcessNum = min(processNum+np, numProcesses)  
            for proc in processes[processNum:endProcessNum]:  
                proc.start()  
            for proc in processes[processNum:endProcessNum]:  
                proc.join()  
            for proc,c in pipe[processNum:endProcessNum]:  
                outputList.append(proc.recv())  

            processNum = endProcessNum  

        threadNum = endThreadNum

    return outputList    



# This is for testing:
def sum2d(arr):
    M, N = arr.shape
    result = 0.0
    for i in range(M):
        for j in range(N):
            result += arr[i,j]
    return result

if __name__ == '__main__':  
    arr = np.ones((1000,1000))
    arr2=[arr+x for x in range(1000)]
    print parmap(sum2d,arr2)
#    print parmapLimitThreads(sum2d,arr2)
#    print parmap(lambda x:x**x,range(1,5))


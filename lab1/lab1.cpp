#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <omp.h>

#define threadsCount 2
static long cntSteps=1000000;

int getTimeExecution(double(*)(double, double), double, double);
int getTimeExecution(double(*op)(double, double), double a, double b){
    int i;
    double cntIteration = 5;
    double totalTime = 0.0;
    for (i=0; i<cntIteration; i++)
    {
        auto t1 = std::chrono::high_resolution_clock::now();
        op(a, b);
        auto t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::micro> fp_ms = t2 - t1;
        totalTime += fp_ms.count();
    }
    return totalTime/cntIteration;
}

double f(double x){
    return pow(sin(1.0/x),2)*pow(1.0/x,2);
}

double integrateSerial(double a, double b){
    double x, temp = 0.0;
    int i;
    double step=(double)(b-a)/cntSteps;

    for (i=1; i<cntSteps; i++)
    {
        x = a+step*i;
        temp += f(x);
    }
    return step * (f(a)/2.0 + temp + f(b)/2.0);
}

double integrateAtomic(double a, double b){
    double x, temp = 0.0;
    int i;
    double step=(double)(b-a)/cntSteps;
    #pragma omp parallel for private(i, x) num_threads(threadsCount)
        for (i=1; i<cntSteps; i++) {
            x = a+step*i;   
            #pragma omp atomic
            temp += f(x);
        }
    return step * (f(a)/2.0 + temp + f(b)/2.0);
}

double integrateCritical(double a, double b){
    double x, temp = 0.0;
    int i;
    double step=(double)(b-a)/cntSteps;
    #pragma omp parallel for private(i, x) num_threads(threadsCount)
        for (i=1; i<cntSteps; i++) {
            x = f(a+step*i);   
            #pragma omp critical
            temp += x;
        }
    return step * (f(a)/2.0 + temp + f(b)/2.0);
}

double integrateLock(double a, double b){
    double x, temp = 0.0;
    int i;
    double step=(double)(b-a)/cntSteps;
    omp_lock_t mylock;
    omp_init_lock(&mylock);
    #pragma omp parallel for private(i, x) num_threads(threadsCount)
        for (i=1; i<cntSteps; i++) {
            x = a+step*i;   
            omp_set_lock(&mylock);
            temp += f(x);
            omp_unset_lock(&mylock);
        }
    omp_destroy_lock(&mylock);
    return step * (f(a)/2.0 + temp + f(b)/2.0);
}

double integrateReduction(double a, double b){
    double x, temp = 2.0;
    int i;
    double step=(double)(b-a)/cntSteps;
        #pragma opm parallel reduction (+ : temp)
        for (i=1; i<cntSteps; i++) {
            temp += f(a+step*i);
        }
    return step * (f(a)/2.0 + temp + f(b)/2.0);
}

int main() {
    double a, b, i = 0.0;
    for (i = 0.00001; i<11; i=i*10) 
    {
        a = i;
        b = i*10;
        double integralValueSerial = integrateSerial(a,b);
        double timeExecutionSerial = getTimeExecution(integrateSerial, a, b);

        double integralValueAtomic = integrateAtomic(a,b);
        double timeExecutionAtomic = getTimeExecution(integrateAtomic, a, b);

        double integralValueCritical = integrateCritical(a,b);
        double timeExecutionCritical = getTimeExecution(integrateCritical, a, b);

        double integralValueLock = integrateLock(a,b);
        double timeExecutionLock = getTimeExecution(integrateLock, a, b);

        double integralValueReduction = integrateReduction(a,b);
        double timeExecutionReduction = getTimeExecution(integrateReduction, a, b);

        double integralEntirely = (2.0*(b-a)/(a*b) + sin(2.0/b) - sin(2.0/a)) / 4.0;
        double accurancy = abs(integralEntirely - integralValueSerial)*100.0/integralEntirely;

        std::cout<<std::endl<< "A \t| B \t| Npoint|" << std::endl;
        std::cout<<" "<<a<<" \t|" <<b<<" \t| "<<cntSteps<<"\t|"<< std::endl;
        
        std::cout<<"----------------|-------|"<<std::endl;
        std::cout<< "Method \t\t|  Value \t| AvrTime  |" << std::endl;
        std::cout<< "Serial\t\t| "<<std::setw(10)<< std::setprecision(11) << integralValueSerial << "|"
                 <<std::setw(10) << std::setprecision(10) << timeExecutionSerial <<"| "<<std::endl;
        std::cout<< "Atomic\t\t| "<<std::setw(10)<< std::setprecision(10) << integralValueAtomic << "|"
                 <<std::setw(10) << std::setprecision(10) << timeExecutionAtomic <<"| "<< std::endl;
        std::cout<< "Critical\t| "<<std::setw(10)<< std::setprecision(10) << integralValueCritical << "|"
                 <<std::setw(10) << std::setprecision(10) << timeExecutionCritical <<"| "<< std::endl;
        std::cout<< "Lock\t\t| "<<std::setw(10)<< std::setprecision(10) << integralValueLock << "|"
                 <<std::setw(10) << std::setprecision(10) << timeExecutionLock <<"| "<< std::endl;
        std::cout<< "Reduction\t| "<<std::setw(10)<< std::setprecision(10) << integralValueReduction << "|"
                 <<std::setw(10) << std::setprecision(10) << timeExecutionReduction <<"| "<< std::endl;
        std::cout << "W/o steps \t| " << std::setw(10) << std::setprecision(10) << integralEntirely << "|"<<std::endl;
        std::cout << "Accurancy \t| "  << std::setw(10) << std::setprecision(11) << accurancy<< "\t |" << std::endl;
    }
}
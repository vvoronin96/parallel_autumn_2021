#include<iostream>
#include<stdlib.h>
#include<time.h>
#include <chrono>
#include <stdio.h>
using namespace std;
static int n=415;

int main(int argc, char **argv)
{
    srand(time(0)); //set state of generate depends on current time

    int firstMatrix[n][n], secondMatrix[n][n], mult[n][n];
    int r1, c1, r2, c2, i, j, k;
    

    // If column of first matrix in not equal to row of second matrix,
    // ask the user to enter the size of matrix again.
    /*while (c1!=r2)
    {
        cout << "Error! column of first matrix not equal to row of second.";

        cout << "Enter rows and columns for first matrix: ";
        cin >> r1 >> c1;

        cout << "Enter rows and columns for second matrix: ";
        cin >> r2 >> c2;
    }*/

    for (int w=10; w<415; w = w+50)
    {
        
        r1=w;
        c1 = w; 
        r2 = c1;
        c2=w;
        // Storing elements of first matrix.
        for(i = 0; i < r1; ++i){
            for(j = 0; j < c1; ++j)
            {
                firstMatrix[i][j] = rand() %10;
            }
        }

        // Storing elements of second matrix.
        for(i = 0; i < r2; ++i){
            for(j = 0; j < c2; ++j)
            {
                secondMatrix[i][j] = rand() %10;
            }
        }
        // Initializing elements of matrix mult to 0.
        for(i = 0; i < r1; ++i)
            for(j = 0; j < c2; ++j)
            {
                mult[i][j]=0;
            }

        ////////////////////Serial multiply///////////////////

        // Measurement time spent and multiplying matrix a and b and storing in array mult.
        double cntIteration = 5;
        double totalTime = 0.0;
        for (int q=0; q<cntIteration; q++)
        {
            auto t1 = std::chrono::high_resolution_clock::now();
            for(i = 0; i < r1; ++i)
            for(j = 0; j < c2; ++j)
                for(k = 0; k < c1; ++k)
                {
                    mult[i][j] += firstMatrix[i][k] * secondMatrix[k][j];
                }
            auto t2 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::micro> fp_ms = t2 - t1;
            totalTime += fp_ms.count();
        }
        double timeSpentSerial = totalTime/cntIteration;
        
        // Displaying the multiplication of two matrix.
        cout << w<<"\t" << "timeSpentSerial: " << timeSpentSerial<<endl;
        /*cout << endl << "Output Matrix: " << endl;
        for(i = 0; i < r1; ++i){
            for(j = 0; j < c2; ++j)
            {
                cout << " " << mult[i][j];
            }
            cout << endl;
        }*/
    }
    return 0;
}
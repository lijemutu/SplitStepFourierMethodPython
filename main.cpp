#include <iostream>
#include <vector>

int main(){


    //int numValues = 128;
    double vecStart = -1;
    double vecFinish = 1;
    double vecInterval = 0.05;
    int numValues = (vecFinish-vecStart)/vecInterval + vecInterval;
    std::vector <double> mainVector(numValues+1);

    for (int i = 0; i < mainVector.size(); i++)
    {
        mainVector[i] = vecStart+vecInterval*i;
        std::cout << mainVector[i] << " ";
        
    }

    return 0; 
}






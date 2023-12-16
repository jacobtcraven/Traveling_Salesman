#include "classes.h"
using namespace std;


int main() {

    string graphName;

    for(int i = 1; i < 11; i++){

        graphName = "graph" + to_string(i) + ".txt";
        bruteForceTSP(graphName);
    
    }
    return 0;
}
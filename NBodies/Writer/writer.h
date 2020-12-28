//
// Created by Максим on 28.12.2020.
//

#ifndef NBODIES_WRITER_H
#define NBODIES_WRITER_H

#include <string>
#include <vector>
#include <fstream>
#include <iomanip>

template<typename Type>
class writer {
public:
    void writeRes(std::string filename, std::vector<Type>& x, std::vector<Type> &y);
};

template <typename Type>
void writer<Type>::writeRes(std::string filename, std::vector<Type> &x, std::vector<Type> &y) {
    std::ofstream newFile;
    newFile.open(filename);
    for(int i = 0; i<x.size(); i++){
        newFile << std::setprecision(64);
        newFile << x[i] << "," << y[i] << std::endl;
    }
    newFile.close();
}


#endif //NBODIES_WRITER_H

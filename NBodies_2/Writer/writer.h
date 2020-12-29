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
    void writeBodies(std::string filename, std::vector<Type> &x_1, std::vector<Type> &y_1, std::vector<Type> &z_1,
                     std::vector<Type> &x_2, std::vector<Type> &y_2, std::vector<Type> &z_2);
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

template <typename Type>
void writer<Type>::writeBodies(std::string filename, std::vector<Type> &x_1, std::vector<Type> &y_1, std::vector<Type> &z_1,
                 std::vector<Type> &x_2, std::vector<Type> &y_2, std::vector<Type> &z_2){
    std::ofstream newFile;
    newFile.open(filename);
    for(int i = 0; i<x_1.size(); i++){
        newFile << std::setprecision(64);
        newFile << x_1[i] << "," << y_1[i] << "," << z_1[i] << "," << x_2[i] << "," << y_2[i] << "," << z_2[i] <<  std::endl;
    }
    newFile.close();
}


#endif //NBODIES_WRITER_H

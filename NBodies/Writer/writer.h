//
// Created by Максим on 28.12.2020.
//

#ifndef NBODIES_WRITER_H
#define NBODIES_WRITER_H

#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include "json.hpp"

using nlohmann::json;

template<typename Type>
class writer {
public:
    void writeRes(std::string filename, const json &data);
};

template <typename Type>
void writer<Type>::writeRes(std::string filename, const json &data) {
    std::ofstream newFile;
    newFile.open(filename);
    newFile << data;
    newFile.close();
}



#endif //NBODIES_WRITER_H

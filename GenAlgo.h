#include "Matrix.h"
#include <experimental/filesystem>
#include <thread>
#include <algorithm>


#ifndef KR_1_MAG_GENALGO_H
#define KR_1_MAG_GENALGO_H


class GenAlgo {
public:
    static const int SRC_CNT = 17;
    static const int MAX_CNT = 100;
private:
    vector<Matrix> data;
    vector<double> values;
    vector<string> filenames;
    vector<pair<double, Matrix>> genitor_data;
    int min_value;
    int max_value;
private:

public:
    GenAlgo()= default;
public:
    void read(const string&);
    void rationing();
    void choose_files();
    double calc_matrix_value(int num);
    void make_100_matrix();
    void calc_source_values();
    void read_filenames(const string&);
    void read_source_values(const string&);
    void make_min_max();
    void Genitor();
public:
    friend std::ostream& operator<< (std::ostream &, const GenAlgo &);
    void write(string);
};


#endif //KR_1_MAG_GENALGO_H

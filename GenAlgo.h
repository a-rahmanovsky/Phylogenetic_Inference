#include "Matrix.h"
#include <experimental/filesystem>
#include <thread>
#include <algorithm>


#ifndef KR_1_MAG_GENALGO_H
#define KR_1_MAG_GENALGO_H


class GenAlgo {
public:
    static const int SRC_CNT = 17;
    int MAX_CNT = 200;
private:
    vector<Matrix> data;
    vector<double> values;
    vector<string> filenames;
    vector<pair<double, Matrix>> genitor_data;
    vector<int> random_data;
    int min_value{};
    int max_value{};
private:
    bool find_matrix(const Matrix&);
    void assign_random_data();
    int get_random_matrix();
public:
    GenAlgo()= default;
    explicit GenAlgo(int max_cnt): MAX_CNT(max_cnt) {};
public:
    void read(const string&, int);
    void rationing();
    void choose_files();
    double calc_matrix_value(int num);
    void make_all_matrix();
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

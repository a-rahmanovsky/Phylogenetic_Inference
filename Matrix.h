#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

#ifndef KR_1_MAG_MATRIX_H
#define KR_1_MAG_MATRIX_H


class Matrix{
private:
    static const int SIZE = 20;
    constexpr static const double RATION_COEF = 4.0;
private:
    double matrix[SIZE][SIZE]{};
    int stable = 10;
    string fields = "CFYWMLIVGPATSNHQEDRK";
public:
    Matrix()= default;
    explicit Matrix (std::string &);
public:
    void rationing();
    void mutation(int min, int rigth);
    double get_min();
    double get_max();
public:
    static Matrix hybridization(const Matrix&, const Matrix&);
public:
    int operator () (int, int) const;
    friend std::ostream& operator<< (std::ostream &, const Matrix &);
    friend std::istream & operator>> (std::istream &, Matrix &);
    Matrix& operator = (const Matrix&);
};


#endif //KR_1_MAG_MATRIX_H

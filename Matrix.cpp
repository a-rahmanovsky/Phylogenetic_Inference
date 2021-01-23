#include "Matrix.h"

std::ostream &operator<<(std::ostream &out, const Matrix &matrix1) {
    for (int i = 0; i < Matrix::SIZE; i++)
        out << '\t' << matrix1.fields[i];
    out << endl;
    for (int i = 0; i < Matrix::SIZE; i++){
        out << matrix1.fields[i] << '\t';
        for (int j = 0; j < Matrix::SIZE; j++){
            out << matrix1.matrix[i][j] << '\t';
        }
        out << std::endl;
    }
    return out;
}

int Matrix::operator()(int i, int j) const {
    return matrix[i][j];
}

std::istream &operator>>(std::istream &in, Matrix &matrix1) {
    for (int i = 0; i < Matrix::SIZE; i++){
        for (int j = 0; j < Matrix::SIZE; j++)
            in >> matrix1.matrix[i][j];
    }
    return in;
}

Matrix::Matrix(std::string &path) {
    std::fstream f{path, ios::in};
    char x;
    for (int i = 0; i < Matrix::SIZE; i++) {
        f >> x;
    }
    for (int i = 0; i < Matrix::SIZE; i++){
        f >> x;
        for (int j = 0; j < Matrix::SIZE; j++) {
            f >> matrix[i][j];
        }
    }
}

void Matrix::rationing() {
    for (int i = 0; i < Matrix::SIZE; i++){
        for (int j = 0; j < Matrix::SIZE; j++){
             if (i == stable && j == stable)
                continue;
            matrix[i][j] = round(matrix[i][j] * (Matrix::RATION_COEF / matrix[stable][stable]));
        }
    }
    matrix[stable][stable] = 4;
}

Matrix Matrix::hybridization(const Matrix &a, const Matrix &b) {
    Matrix newm;
    for (int i = 0; i < Matrix::SIZE; i++){
        for (int j = i; j < Matrix::SIZE; j++){
            int maxv = max(a.matrix[i][j], b.matrix[i][j]);
            int minv = min(a.matrix[i][j], b.matrix[i][j]);
            int new_value = minv + rand() % (maxv - minv + 1);
            newm.matrix[i][j] = new_value;
            newm.matrix[j][i] = new_value;
        }
    }
    return newm;
}

void Matrix::mutation(int min, int max) {
    int x = rand() % 20;
    for (int y = 0; y < 20; y++){
        if (x == stable && y == stable){
            continue;
        }
        int type = rand() % 2;
        if (type) {
            if (matrix[x][y] < max) {
                matrix[x][y]++;
                matrix[y][x]++;
            }
            else{
                matrix[x][y]--;
                matrix[y][x]--;
            }
        }
        else {
            if (matrix[x][y] > min) {
                matrix[x][y]--;
                matrix[y][x]--;
            }
            else{
                matrix[x][y]++;
                matrix[y][x]++;
            }
        }
    }
}

Matrix &Matrix::operator=(const Matrix &matrix1){
    if (this == &matrix1){
        return *this;
    }
    for (int i = 0; i < SIZE; i++){
        for (int j = 0; j < SIZE; j++){
            matrix[i][j] = matrix1.matrix[i][j];
        }
    }
    return *this;
}

double Matrix::get_min() {
    double min = matrix[0][0];
    for (int i = 0; i < Matrix::SIZE; i++){
        for (int j = 0; j < Matrix::SIZE; j++){
            if (i == stable && j == stable) continue;
            min = fmin(min, matrix[i][j]);
        }
    }
    return min;
}

double Matrix::get_max() {
    double max = matrix[0][0];
    for (int i = 0; i < Matrix::SIZE; i++){
        for (int j = 0; j < Matrix::SIZE; j++){
            if (i == stable && j == stable) continue;
            max = fmax(max, matrix[i][j]);
        }
    }
    return max;
}


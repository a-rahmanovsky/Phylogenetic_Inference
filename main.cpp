#include "GenAlgo.h"


int main()
{
    srand(time(NULL));
    GenAlgo ga;

    // Нормируем исходные
//    ga.read("source_matrix");
//    ga.rationing();
//    ga.write("result_matrix");

    // 17 в 100
//    ga.read("result_matrix");
//    ga.make_100_matrix();
//    ga.write("result_matrix");

    // Считаем исходные значения
//    ga.read("result_matrix");
//    ga.read_filenames("type_data/train_data.txt");
//    ga.calc_source_values();

    ga.read("result_matrix");
    ga.read_filenames("type_data/train_data.txt");
    ga.read_source_values("source_values/all_values.txt");
    ga.make_min_max();
    ga.Genitor();
}
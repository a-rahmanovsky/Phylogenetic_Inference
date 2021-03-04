#include "GenAlgo.h"

void choose_data(GenAlgo &ga){
    ga.choose_files();
}

void normalize_source(GenAlgo& ga){
    ga.read("source_matrix", ga.SRC_CNT);
    ga.rationing();
    ga.write("result_matrix");
}

void create_all_from_17(GenAlgo &ga, int max){
    ga.read("result_matrix", ga.SRC_CNT);
    ga.make_all_matrix();
    ga.write("result_matrix");
}

void eval_base_values(GenAlgo &ga){
    ga.read("result_matrix", ga.MAX_CNT);
    ga.read_filenames("type_data/train_data.txt");
    ga.calc_source_values();
}

void Genitor(GenAlgo &ga){
    ga.read("result_matrix", ga.MAX_CNT);
    ga.read_filenames("type_data/train_data.txt");
    ga.read_source_values("source_values/all_values.txt");
    ga.make_min_max();
    ga.Genitor();
}

int main()
{
    srand(time(NULL));
    GenAlgo ga(300);

    // Выбираем train/test/validate
    //choose_data(ga);

   //  Нормируем исходные
    //normalize_source(ga);

    // 17 в all
    //create_all_from_17(ga, 300);

    // Считаем исходные значения
    //eval_base_values(ga);

    // Запускаем Genitor
     Genitor(ga);

}
//
// Created by andrey on 24.12.2020.
//

#include "GenAlgo.h"

#include <random>

std::ostream &operator<<(ostream &out, const GenAlgo &algo) {
    for (int i = 1; i <= GenAlgo::SRC_CNT; i++){
        out << "number: " << i << endl;
        out << algo.data[i - 1];
        out << endl;
    }
    return out;
}

void GenAlgo::read(const string& source_path, int cnt) {
    for (int i = 1; i <= cnt; i++){
        string path = source_path + "/" + to_string(i) + ".txt";
        Matrix matrix(path);
        data.push_back(matrix);
    }
}

void GenAlgo::rationing() {
    for (int i = 0; i < data.size(); i++){
        data[i].rationing();
    }
}

void GenAlgo::write(string path) {
    ofstream f;
    for (int i = 0; i < data.size(); i++){
        f.open(path + '/' + to_string(i + 1) + ".txt");
        f << data[i];
        f.close();
    }
}

void GenAlgo::choose_files() {
    std::string path = "data/Alignments";
    vector<string> all_files;
    for (auto & p : std::experimental::filesystem::directory_iterator(path)) {
        string name = "";
        for (int i = path.size() + 1; i < p.path().string().size() - 4; i++)
            name += p.path().string()[i];
        all_files.push_back(name);
    }
    cout << all_files.size() << endl;
    fstream f;
    f.open("type_data/train_data.txt", ios::out | ios::ate);
    for (int i = 0; !all_files.empty(); i++){
        int pos = rand() % all_files.size();
        f << all_files[pos] << endl;
        all_files.erase(all_files.begin() + pos);
    }
    f.close();
    f.open("type_data/validate_data.txt", ios::out | ios::ate);
    for (int i = 0; i < 300; i++){
        int pos = rand() % all_files.size();
        f << all_files[pos] << endl;
        all_files.erase(all_files.begin() + pos);
    }
    f.close();
    f.open("type_data/test_data.txt", ios::out | ios::ate);
    for (int i = 0; i < all_files.size(); i++){
        f << all_files[i] << endl;
    }
    f.close();
}

void exec_thread(int begin, int end, vector<string> &filenames, string id, string num_m) {
    for (int i = begin; i <= end; i++) {
        string name = filenames[i];
        string c1 = "pq/pq -alignment data/Alignments/" + name;
        c1 += ".afa -pwm result_matrix/" + num_m + ".txt -out trees/res" + id + ".txt";
        c1 += " -grType one -nniType none -randLeaves 0";

        string c2 = "rf_dist/a.out data/References/" + name;
        c2 += ".tre trees/res" + id + ".txt";
        system(c1.c_str());
        system(c2.c_str());
    }
}

double GenAlgo::calc_matrix_value(int num) {
    auto begin = std::chrono::steady_clock::now();
    int cnt_threads = 10;
    int size = filenames.size() / cnt_threads;
    vector<thread> threads;
    for (int i = 0; i < cnt_threads; i++) {
        threads.emplace_back(exec_thread, size * i, size * (i + 1) - 1,
               ref(filenames), to_string(i + 1),
               to_string(num));
    }
    for (int i = 0; i < threads.size(); i++){
        threads[i].join();
    }
    fstream f{"result.txt"};
    double sum = 0;
    while (!f.eof()){
        double value;
        f >> value;
        sum += value;
    }
    auto end = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::seconds>(end - begin);
    cout << "time: " << elapsed_ms.count() << endl;
    ofstream {"result.txt"};
    return sum / filenames.size();
}

void GenAlgo::make_all_matrix() {
    for (int i = 18; i <= MAX_CNT; i++){
        int first = 0;
        int second = rand() % 18;
        while (first == second) second = rand() % 18;
        Matrix newm = Matrix::hybridization(data[first], data[second]);
        newm.mutation(min_value, max_value);
        data.push_back(newm);
    }
}

void GenAlgo::calc_source_values() {
    ofstream f{"source_values/all_values.txt", ios::out | ios::ate};
    for (int i = 1; i <= MAX_CNT; i++){
        double res = calc_matrix_value(i);
        f << res << endl;
        cout << "number: " << i << ", value: " << res << endl;
    }
}

void GenAlgo::read_filenames(const string& path) {
    fstream f{path, ios::in};
    while (!f.eof()){
        string name;
        f >> name;
        filenames.push_back(name);
    }
    filenames.pop_back();
}

bool comp(pair<double, Matrix> a, pair<double, Matrix> b){
    return a.first < b.first;
}

void GenAlgo::read_source_values(const string& path) {
    fstream f{path, ios::in};
    int pos = 0;
    while (pos < MAX_CNT){
        double value;
        f >> value;
        genitor_data.emplace_back(value, data[pos]);
        pos++;
    }
    sort(genitor_data.begin(), genitor_data.end(), comp);
}

void GenAlgo::Genitor() {
    assign_random_data();
    int cnt = 0;
    int type = 1;
    Matrix newm;
    cout << data.size() << endl;
    while (true){
        if (cnt % 20 == 0){
            write("backup_genitor_matrix");
            cout << "Current best value: " << genitor_data.front().first << endl;
        }
        if (type == 1){
            int first_pos = get_random_matrix(), second_pos = get_random_matrix();
            auto first = genitor_data[first_pos].second;
            auto second = genitor_data[second_pos].second;
            newm = Matrix::hybridization(first, second);
            cout << "Hybridization: " << first_pos << " and " << second_pos << endl;
        }
        else{
            int pos = get_random_matrix();
            auto source_matrix = genitor_data[pos].second;
            newm = source_matrix;
            newm.mutation(min_value, max_value);
            cout << "Mutation: " << pos << endl;
        }
        ofstream f{"result_matrix/-1.txt"};
        f << newm;
        double res = calc_matrix_value(-1);
        cout << "Value: " << res << ", worst value: " << genitor_data.back().first << endl;
        pair<double, Matrix> new_element = {res, newm};
        if (res < genitor_data.back().first && !find_matrix(newm)){
            genitor_data.pop_back();
            auto pos = genitor_data.begin();
            for (; pos->first < res; pos++);
            genitor_data.insert(pos, new_element);
            cout << "Insert in " << pos - genitor_data.begin() << " position" << endl;
            if (pos == genitor_data.begin()){
                cnt = 0;
                cout << "New best matrix with res: " << res << endl;
                ofstream fbest{"best_matrix.txt"};
                fbest << newm;
                fbest << endl;
                fbest << res;
            }
            else {
                cnt++;
                cout << "Current cnt: " << cnt << endl;
            }
        }
        else{
            cout << "Bad matrix" << endl;
        }
        type *= -1;
    }
}

void GenAlgo::make_min_max() {
    double min = data[0].get_min();
    double max = data[0].get_max();
    for (int i = 1; i < data.size(); i++){
        min = fmin(min, data[i].get_min());
        max = fmax(max, data[i].get_max());
    }
    max_value = int(max);
    min_value = int(min);
}

bool GenAlgo::find_matrix(const Matrix& matrix) {
    for (auto element : genitor_data){
        if (element.second == matrix)
            return true;
    }
    return false;
}

void GenAlgo::assign_random_data() {
    for (int i = 0; i < MAX_CNT; i++){
        for (int j = 0; j < MAX_CNT - i; j++)
            random_data.push_back(i);
    }
    shuffle(random_data.begin(), random_data.end(), std::mt19937(std::random_device()()));
}

int GenAlgo::get_random_matrix() {
    int pos = rand() % random_data.size();
    //cout << random_data.size() << " " << pos << endl;
    return random_data[pos];
}





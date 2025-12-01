#include "counting.h"
#include <cstdlib>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <epsilon> <dataset> [num_rounds] [vertex_ratio] [algo_switch]" << std::endl;
        std::cout << "Example: " << argv[0] << " 1.0 graphs/email.edges 40 0.1 0" << std::endl;
        return 1;
    }

    double eps = std::stod(argv[1]);
    std::string dataset = argv[2];
    int rounds = (argc > 3) ? std::atoi(argv[3]) : 40;
    double v_ratio = (argc > 4) ? std::stod(argv[4]) : 0.1;
    int algo_sw = (argc > 5) ? std::atoi(argv[5]) : 0;

    run_kcore_suite(eps, dataset, rounds, v_ratio, algo_sw);
    return 0;
}

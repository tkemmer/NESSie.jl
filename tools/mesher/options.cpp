#include "options.h"
#include <iostream>

void printUsage() {
    std::cout << "\n\033[1mProteinES.jl molecule mesher (using GAMer)\033[0m\n"
         << "==========================================\n\n"
         << "\033[1mUsage:\033[0m\n\n"
         << "\tmesher <file in> <dir out> [<option>...]\n"
         << "\n\033[1mSupported input formats:\033[0m\n\n"
         << "\t\033[1mPDB\033[0m\t*.pdb\n"
         << "\t\033[1mPQR\033[0m\t*.pqr\n"
         << "\n\033[1mSupported options (with default values):\033[0m\n\n";
    for(const auto& e: options)
        std::cout << "\t\033[1m--" << e.first << "=<" << e.second.type << ">\033[0m\t"
             << e.second.defval << "\t" << e.second.desc << "\n";
    std::cout << std::endl;
}

Options parseOptions(char** first, char** last) {
    Options opt;
    for(auto it = first; it != last; ++it) {
        std::string curr {(*it)+2};   // omit "--"
        auto p = curr.find("="); // split param and value
        // find corresponding option and call the init function
        options.at(curr.substr(0, p)).init(opt, (*it)+p+3);
    }
    return opt;
}

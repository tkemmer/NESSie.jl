#include "options.h"
#include <iostream>

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

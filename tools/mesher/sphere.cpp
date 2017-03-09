#include <iostream>
#include "functions.h"

using namespace std;

void printUsage() {
    std::cout << "\n\033[1mProteinES.jl sphere generator (using GAMer)\033[0m\n"
        << "===========================================\n\n"
        << "Generates a surface mesh for a sphere with the given radius.\n\n"
        << "\033[1mUsage:\033[0m\n\n"
        << "\tsphere <file out> <sphere radius> <sphere quality>\n"
        << "\n\033[1mSphere quality:\033[0m\n\n"
        << "\tDetermines the number of vertices (4^n + 2)\n\n"
        << std::endl;
}

int main(int argc, char** argv) {
    if(argc != 4) {
        printUsage();
        return 1;
    }

    Options opt;
    opt.sphere_ratio   = std::atof(argv[2]);
    opt.sphere_quality = std::atoi(argv[3]);

    cout << "\033[1mGenerating sphere...\033[0m\n";
    auto surf = generateBoundingSphere({0, 0, 0, 1}, opt);

    cout << "\033[1mWriting surface mesh...\033[0m\n";
    SurfaceMesh_writeOFF(surf.data(), argv[1]);
}

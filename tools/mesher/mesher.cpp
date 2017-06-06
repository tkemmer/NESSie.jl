#include <iostream>
#include "functions.h"

using namespace std;

void printUsage() {
    cout << "\n\e[1mNESSie.jl molecule mesher (using GAMer)\e[0m\n"
         << "==========================================\n\n"
         << "\e[1mUsage:\e[0m\n\n"
         << "\tmesher <file in> <dir out> [<option>...]\n"
         << "\n\e[1mSupported input formats:\e[0m\n\n"
         << "\t\e[1mPDB\e[0m\t*.pdb\n"
         << "\t\e[1mPQR\e[0m\t*.pqr\n"
         << "\n\e[1mSupported options (with default values):\e[0m\n\n";
    for(const auto& e: options)
        cout << "\t\e[1m--" << e.first << "=<" << e.second.type << ">\e[0m\t"
             << e.second.defval << "\t" << e.second.desc << "\n";
    cout << endl;
}

int main(int argc, char** argv) {
    if(argc < 3) {
        printUsage();
        return 1;
    }
    string pdb_in {argv[1]};
    string out_dir {argv[2]};
    auto opt = parseOptions(&argv[3], &argv[argc]);

    cout << "\e[1mGenerating surface mesh...\e[0m\n";
    Surface molsurf;
    AtomList atoms;
    tie(molsurf, atoms) = gaussianSurfaceFromPDB(pdb_in, opt.iso_value);

    cout << "\e[1mSmoothen surface mesh...\e[0m\n";
    smoothSurfaceMesh(molsurf, opt);

    cout << "\e[1mWriting surface mesh...\e[0m\n";
    SurfaceMesh_writeOFF(molsurf.data(), &(out_dir + "/mesh.off")[0]);

    cout << "\e[1mGenerating bounding sphere...\e[0m\n";
    auto center = SurfaceMesh_getCenterRadius(molsurf.data());
    auto sphere = generateBoundingSphere(center, opt);
    cout << "Molecule radius: " << center.radius << ", Center: ["
         << center.x << ", " << center.y << ", " << center.z << "]\n";

    cout << "\e[1mGenerating volume mesh...\e[0m\n";
    Surface mesh = SurfaceMesh_merge(molsurf.data(), sphere.data());
    auto domain  = generateVolumeMesh(mesh, atoms);

    cout << "\e[1mWriting volume mesh...\e[0m\n";
    cout << "\e[1m * Inner mesh...\e[0m\n";
    writeMcsf(out_dir, *domain, true, opt.sphere_ratio, center);
    cout << "\e[1m * Outer mesh...\e[0m\n";
    writeMcsf(out_dir, *domain, false, opt.sphere_ratio, center);
}

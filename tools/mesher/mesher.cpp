#include <iostream>
#include "functions.h"

using namespace std;

void printUsage() {
    cout << "\n\033[1mProteinES.jl molecule mesher (using GAMer)\033[0m\n"
         << "==========================================\n\n"
         << "\033[1mUsage:\033[0m\n\n"
         << "\tmesher <file in> <dir out> [<option>...]\n"
         << "\n\033[1mSupported input formats:\033[0m\n\n"
         << "\t\033[1mPDB\033[0m\t*.pdb\n"
         << "\t\033[1mPQR\033[0m\t*.pqr\n"
         << "\n\033[1mSupported options (with default values):\033[0m\n\n";
    for(const auto& e: options)
        cout << "\t\033[1m--" << e.first << "=<" << e.second.type << ">\033[0m\t"
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

    cout << "\033[1mGenerating surface mesh...\033[0m\n";
    Surface molsurf;
    AtomList atoms;
    tie(molsurf, atoms) = gaussianSurfaceFromPDB(pdb_in, opt.iso_value);

    cout << "\033[1mSmoothen surface mesh...\033[0m\n";
    smoothSurfaceMesh(molsurf, opt);

    cout << "\033[1mWriting surface mesh...\033[0m\n";
    SurfaceMesh_writeOFF(molsurf.data(), &(out_dir + "/mesh.off")[0]);

    cout << "\033[1mGenerating bounding sphere...\033[0m\n";
    auto center = SurfaceMesh_getCenterRadius(molsurf.data());
    auto sphere = generateBoundingSphere(center, opt);
    cout << "Molecule radius: " << center.radius << ", Center: ["
         << center.x << ", " << center.y << ", " << center.z << "]\n";

    cout << "\033[1mGenerating volume mesh...\033[0m\n";
    Surface mesh = SurfaceMesh_merge(molsurf.data(), sphere.data());
    auto domain = generateVolumeMesh(mesh, atoms);

    cout << "\033[1mWriting volume mesh...\033[0m\n";
    cout << "\033[1m * Inner mesh...\033[0m\n";
    writeMcsf(out_dir, *domain, true, opt.sphere_ratio, center);
    cout << "\033[1m * Outer mesh...\033[0m\n";
    writeMcsf(out_dir, *domain, false, opt.sphere_ratio, center);
}

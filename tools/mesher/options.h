#ifndef MESHER_OPTIONS_H
#define MESHER_OPTIONS_H

#include <cstdlib>
#include <functional>
#include <map>
#include <string>

struct Options {
    int sphere_quality {6};     // determines the number of vertices of the bounding sphere (4^n + 2)
    float sphere_ratio {2.0f};  // sphere-to-molecule size ratio
    int max_iterations {6};     // maximum number of iterations for surface smoothing
    int mesh_nodes_min {10000}; // minimum number of mesh vertices
    int mesh_nodes_max {80000}; // maximum number of mesh vertices
    float coarsen_rate {0.3f};  // Coarsening rate
    float iso_value    {1.4f};  // used for gaussian surface
};

struct Option {
    std::string type;
    std::string defval;
    std::string desc;
    std::function<void(Options&, char*)> init;
};

const std::map<std::string, Option> options {
    {"sphere-quality", {
        "int", "6",
        "Determines the number of bounding sphere vertices (4^n + 2)",
        [](Options& opt, char* val) { opt.sphere_quality = std::atoi(val); }
    }},
    {"sphere-ratio", {
        "float", "2.0",
        "Sphere-to-molecule size ratio",
        [](Options& opt, char* val) { opt.sphere_ratio = std::atof(val); }
    }},
    {"max-iterations", {
        "int", "6",
        "Maximum number of surface smoothing steps",
        [](Options& opt, char* val) { opt.max_iterations = std::atoi(val); }
    }},
    {"mesh-nodes-min", {
        "int", "10000",
        "Minimum number of surface mesh vertices",
        [](Options& opt, char* val) { opt.mesh_nodes_min = std::atoi(val); }
    }},
    {"mesh-nodes-max", {
        "int", "80000",
        "Maximum number of surface mesh vertices",
        [](Options& opt, char* val) { opt.mesh_nodes_max = std::atoi(val); }
    }},
    {"coarsen-rate", {
        "float", "0.3",
        "Coarsening rate",
        [](Options& opt, char* val) { opt.coarsen_rate = std::atof(val); }
    }},
    {"iso-value", {
        "float", "1.4",
        "Iso value for gaussian surface (marching cube)",
        [](Options& opt, char* val) { opt.iso_value = std::atof(val); }
    }}
};

void printUsage();
Options parseOptions(char** first, char** last);

#endif // MESHER_OPTIONS_H

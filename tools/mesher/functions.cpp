#include "functions.h"
#include <iostream>

using namespace std;

tuple<Surface, AtomList> gaussianSurfaceFromPDB(string& pdb_in, float iso_value) {
    cout << "\033[1m * Blurring PDB coordinates...\033[0m\n";
    AtomList atoms;
    float* dataset {nullptr};
    int xdim {0}, ydim {0}, zdim {0};
    float min[3], max[3];
    auto max_density = PDB2Volume(&pdb_in[0], &dataset, &xdim, &ydim, &zdim, min, max, &atoms.data, &atoms.size, false);
    float span[] {
        (max[0] - min[0]) / (xdim - 1),
        (max[1] - min[1]) / (ydim - 1),
        (max[2] - min[2]) / (zdim - 1)
    };

    cout << "\033[1m * Extracting isosurfaces...\033[0m\n";
    iso_value = std::min(iso_value, .44f * max_density);
    SamplePoint* holes;
    auto mesh = SurfaceMesh_marchingCube(xdim, ydim, zdim, dataset, iso_value, &holes);
    free(dataset);
    free(holes);

    cout << "\033[1m * Adjusting vertex positions...\033[0m\n";
    // convert pixel -> angstrom
    for(auto i = 0; i < mesh->nv; ++i) {
        auto& vertex = mesh->vertex[i];
        vertex.x = vertex.x * span[0] + min[0];
        vertex.y = vertex.y * span[1] + min[1];
        vertex.z = vertex.z * span[2] + min[2];
    }

    return make_tuple(move(mesh), move(atoms));
}

void smoothSurfaceMesh(Surface& surface, const Options& opt) {
    auto mesh = surface.data();
    SurfaceMesh_createNeighborlist(mesh);

    // smooth surface by moving the surface vertices only
    SurfaceMesh_smooth(mesh, 10, 160, opt.max_iterations, false);

    do {
        // Refine mesh
        if(mesh->nv < opt.mesh_nodes_min) {
            SurfaceMesh_refine(mesh);
            SurfaceMesh_smooth(mesh, 10, 150, opt.max_iterations, true);
            continue;
        }

        // coarsen mesh
        SurfaceMesh_coarse(mesh, opt.coarsen_rate, 1.0, 1.0, 0.5);
        for(auto i = 0; i < opt.max_iterations; ++i) // this seems to preserve features
            if(SurfaceMesh_smooth(mesh, 20, 140, 1, true))
                break;
        SurfaceMesh_normalSmooth(mesh);
    } while(mesh->nv < opt.mesh_nodes_min || mesh->nv > opt.mesh_nodes_max);

    SurfaceMesh_destroyNeighborlist(mesh);
}

Surface generateBoundingSphere(const Atom& center, const Options& opt) {
    // Generate unit sphere mesh
    auto sphere = SurfaceMesh_sphere(opt.sphere_quality);
    // Expand and re-center
    for(auto i = 0; i < sphere->nv; ++i) {
        auto& vertex = sphere->vertex[i];
        vertex.x = vertex.x * center.radius * opt.sphere_ratio + center.x;
        vertex.y = vertex.y * center.radius * opt.sphere_ratio + center.y;
        vertex.z = vertex.z * center.radius * opt.sphere_ratio + center.z;
    }
    return sphere;
}

unique_ptr<Volume> generateVolumeMesh(Surface& mesh, AtomList& atoms) {
    Volume in;
    in.firstnumber = 1;

    // copy all vertices
    auto surface = mesh.data();
    in.numberofpoints = surface->nv;
    in.pointlist = new double[in.numberofpoints * 3];
    for(auto i = 0; i < in.numberofpoints; ++i) {
        in.pointlist[3 * i + 0] = surface->vertex[i].x;
        in.pointlist[3 * i + 1] = surface->vertex[i].y;
        in.pointlist[3 * i + 2] = surface->vertex[i].z;
    }

    // add boundary markers
    in.pointmarkerlist = new int[in.numberofpoints];
    for(auto i = 0; i < in.numberofpoints; ++i)
        in.pointmarkerlist[i] = 1;

    // copy all facets
    in.numberoffacets = surface->nf;
    in.facetlist = new Volume::facet[in.numberoffacets];
    for(auto i = 0; i < in.numberoffacets; ++i) {
        auto facet = &in.facetlist[i];
        facet->holelist = nullptr;
        facet->numberofholes = 0;
        facet->numberofpolygons = 1;
        facet->polygonlist = new Volume::polygon[facet->numberofpolygons];
        auto polygon = &facet->polygonlist[0];
        polygon->numberofvertices = 3;
        polygon->vertexlist = new int[polygon->numberofvertices];
        polygon->vertexlist[0] = surface->face[i].a + in.firstnumber;
        polygon->vertexlist[1] = surface->face[i].b + in.firstnumber;
        polygon->vertexlist[2] = surface->face[i].c + in.firstnumber;
    }

    // insert atoms as nodes
    in.numberofaddpoints = atoms.size;
    in.addpointlist = new double[in.numberofaddpoints * 3];
    for(auto i = 0; i < atoms.size; ++i) {
        in.addpointlist[3 * i + 0] = atoms.data[i].x;
        in.addpointlist[3 * i + 1] = atoms.data[i].y;
        in.addpointlist[3 * i + 2] = atoms.data[i].z;
    }

    // generate Volume mesh
    auto out = new Volume;
    string switches {"npq1.333AAYY"};
    tetrahedralize(&switches[0], &in, out, nullptr);
    return unique_ptr<Volume>(out);
}

void writeMcsf(const string& out_dir, Volume& mesh, bool inner, float sphere_ratio, const Atom& center) {
    auto out = GemMesh_fromPdb(
        &mesh,
        center.radius * sphere_ratio,
        center.x,
        center.y,
        center.y,
        nullptr,
        inner ? 2 : 3
    );
    GemMesh_writeMcsf(out, &(out_dir + "/out." + (inner ? "in" : "out") + ".m")[0]);
    GemMesh_dtor(out);
}

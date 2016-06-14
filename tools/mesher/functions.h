#ifndef MESHER_GAMER_H
#define MESHER_GAMER_H

#include <memory>
#include <tuple>
#include "types.h"
#include "options.h"

std::tuple<Surface, AtomList> gaussianSurfaceFromPDB(std::string& pdb_in, float iso_value);
void smoothSurfaceMesh(Surface& surface, const Options& opt);
Surface generateBoundingSphere(const Atom& center, const Options& opt);
std::unique_ptr<Volume> generateVolumeMesh(Surface& mesh, AtomList& atoms);
void writeMcsf(const std::string& out_dir, Volume& mesh, bool inner, float sphere_ratio, const Atom& center);

#endif // MESHER_GAMER_H

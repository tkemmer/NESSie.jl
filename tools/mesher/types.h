#ifndef MESHER_TYPES_H
#define MESHER_TYPES_H

#include <gamer/gamer.h>
#undef min // Some people just like to watch the world burn.

using Atom = ATOM;
using Volume = tetgenio;

class Surface {
public:
    Surface();
    Surface(SurfaceMesh* mesh);

    Surface(const Surface&) = delete;
    Surface& operator=(const Surface&) = delete;

    Surface(Surface&& surface) noexcept;
    Surface& operator=(Surface&& surface) noexcept;

    ~Surface();

    SurfaceMesh* data() noexcept { return data_; }

private:
    SurfaceMesh* data_;
};

class AtomList {
public:
    AtomList();

    AtomList(const AtomList&) = delete;
    AtomList& operator=(const AtomList&) = delete;

    AtomList(AtomList&& atoms) noexcept;
    AtomList& operator=(AtomList&& atoms) noexcept;

    ~AtomList();

    Atom* data;
    int size;
};

#endif // MESHER_TYPES_H

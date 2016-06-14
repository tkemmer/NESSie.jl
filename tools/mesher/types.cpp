#include "types.h"
#include <utility>

using std::swap;

Surface::Surface():
    data_(nullptr) {
}

Surface::Surface(SurfaceMesh* mesh):
    data_(mesh) {
}

Surface::Surface(Surface&& surface) noexcept:
    data_{nullptr} {
    swap(data_, surface.data_);
}

Surface& Surface::operator=(Surface&& surface) noexcept {
    if(this == &surface)
        return *this;
    swap(data_, surface.data_);
    return *this;
}

Surface::~Surface() {
    if(data_)
        SurfaceMesh_dtor(data_);
}

AtomList::AtomList():
    data(nullptr),
    size(0) {
}

AtomList::AtomList(AtomList&& atoms) noexcept:
    data(nullptr),
    size(atoms.size) {
    swap(data, atoms.data);
}

AtomList& AtomList::operator=(AtomList&& atoms) noexcept {
    if(this == &atoms)
        return *this;
    swap(data, atoms.data);
    swap(size, atoms.size);
    return *this;
}

AtomList::~AtomList() {
    if(data)
        free(data);
}

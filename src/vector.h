#ifndef VECTOR_H
#define VECTOR_H

#include "buffered_reader.h"
#include "buffered_writer.h"
#include "config.h"
#include "point.h"
#include "vector_mixin.h"
#include <array>

template <class R, unsigned D>
struct Vector : VectorMixin<Vector<R, D>, std::array<R, D>, D> {
    using Container = std::array<R, D>;

    Container array{};

    Container &vector_container() { return array; }
    const Container &vector_container() const { return array; }

    R &operator[](unsigned i) { return this->vector_container()[i]; }
    const R &operator[](unsigned i) const
    {
        return this->vector_container()[i];
    }

    unsigned size() const { return D; }
};

template <class R, unsigned D>
void read_varint(BufferedReader &f, Vector<R, D> &vec)
{
    for (auto &x : vec.array) {
        auto v = read_varint(f);
        assert(v >= std::numeric_limits<R>::min());
        assert(v <= std::numeric_limits<R>::max());
        x = static_cast<R>(v);
    }
}

template <class R, unsigned D>
void write_varint(BufferedWriter &f, const Vector<R, D> &vec)
{
    for (const auto &x : vec.array)
        write_varint(f, x);
}

#endif

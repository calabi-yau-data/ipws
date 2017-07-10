#ifndef VECTOR_LIKE_H
#define VECTOR_LIKE_H

#include <iostream>

template <class Base, class Ring, unsigned D>
class VectorLike : public Base {
    using Ret = VectorLike<Base, Ring, D>;

public:
    int compare(const Base &rhs) const {
        for (unsigned i = 0; i < D; ++i) {
            if (Base::vector_like_data(i) < rhs.Base::vector_like_data(i))
                return -1;
            if (Base::vector_like_data(i) > rhs.Base::vector_like_data(i))
                return 1;
        }

        return 0;
    }

    bool operator==(const Base &rhs) const { return compare(rhs) == 0; }
    bool operator!=(const Base &rhs) const { return compare(rhs) != 0; }
    bool operator<(const Base &rhs) const { return compare(rhs) < 0; }
    bool operator>(const Base &rhs) const { return compare(rhs) > 0; }

    Ret operator-() const {
        Ret ret{};
        for (unsigned i = 0; i < D; ++i)
            ret.Base::vector_like_data(i) = -Base::vector_like_data(i);
        return ret;
    }

    Ret operator+(const Base &rhs) const {
        Ret ret{};
        for (unsigned i = 0; i < D; ++i)
            ret.Base::vector_like_data(i) = Base::vector_like_data(i)
                + rhs.Base::vector_like_data(i);
        return ret;
    }

    Ret operator-(const Base &rhs) const {
        Ret ret{};
        for (unsigned i = 0; i < D; ++i)
            ret.Base::vector_like_data(i) = Base::vector_like_data(i)
                - rhs.Base::vector_like_data(i);
        return ret;
    }

    Ret operator*(Ring rhs) const {
        Ret ret{};
        for (unsigned i = 0; i < D; ++i)
            ret.Base::vector_like_data(i) = Base::vector_like_data(i) * rhs;
        return ret;
    }

    friend Ret operator*(Ring lhs, const Base &rhs) {
        Ret ret{};
        for (unsigned i = 0; i < D; ++i)
            ret.Base::vector_like_data(i) = rhs.vector_like_data(i) * lhs;
        return ret;
    }

    Ret operator/(Ring rhs) const {
        Ret ret{};
        for (unsigned i = 0; i < D; ++i)
            ret.Base::vector_like_data(i) = Base::vector_like_data(i) / rhs;
        return ret;
    }

    Ret &operator+=(const Base &rhs) {
        for (unsigned i = 0; i < D; ++i)
            Base::vector_like_data(i) = Base::vector_like_data(i)
                + rhs.Base::vector_like_data(i);
        return *this;
    }

    Ret &operator-=(const Base &rhs) {
        for (unsigned i = 0; i < D; ++i)
            Base::vector_like_data(i) = Base::vector_like_data(i)
                + rhs.Base::vector_like_data(i);
        return *this;
    }

    Ret &operator*=(Ring rhs) {
        for (unsigned i = 0; i < D; ++i)
            Base::vector_like_data(i) *= rhs;
        return *this;
    }

    Ret &operator/=(Ring rhs) {
        for (unsigned i = 0; i < D; ++i)
            Base::vector_like_data(i) /= rhs;
        return *this;
    }

    friend std::ostream &operator<<(std::ostream &os, const Base &rhs) {
        os << "(";
        if (D != 0)
            os << rhs.vector_like_data(0);
        for (unsigned i = 1; i < D; ++i)
            os << ", " << rhs.vector_like_data(i);
        return os << ")";
    }
};

#endif

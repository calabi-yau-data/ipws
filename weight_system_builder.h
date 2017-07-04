#ifndef WEIGHT_SYSTEM_BUILDER
#define WEIGHT_SYSTEM_BUILDER

#include <algorithm>
#include <bitset>
#include <vector>
#include "stl_utils.h"
#include "weight_system.h"
#include "vector.h"

class WeightSystemBuilder {
    struct Generator {
        WeightSystem eq;
        std::bitset<dim> incidences;
    };
    std::vector<Generator> generators; // TODO: can there ever be more generators than dimensions?
    unsigned iteration_;
    // std::array<std::pair<int, int>, dim - 1> symmetries;
    // int num_symmetries;

public:
    // bool operator<(const WeightSystemBuilder &rhs) const
    // {
    //     if (generators.size() != rhs.generators.size())
    //         return generators.size() < rhs.generators.size();

    //     for (size_t i = 0; i < generators.size(); ++i)
    //         if (generators[i].eq != rhs.generators[i].eq)
    //             return generators[i].eq < rhs.generators[i].eq;

    //     return false;
    // }

    WeightSystemBuilder(Long r_numerator, Long r_denominator) : iteration_{0} {
        // num_symmetries = dim - 1;
        // for (int i = 0; i < dim - 1; ++i)
        //     symmetries[i] = {i, i + 1};

        generators.resize(dim);
        for (int generator_nr = 0; generator_nr < dim; ++generator_nr) {
            auto &generator = generators[generator_nr];

            for (int i = 0; i < dim; ++i)
                generator.eq.a[i] = generator_nr == i ? r_numerator : 0;
            generator.eq.c = r_denominator;

            generator.incidences.set(generator_nr);
        }
    }

    unsigned iteration() const { return iteration_; }

    unsigned generator_count() const { return generators.size() ; }

    bool has_symmetry(unsigned idx1, unsigned idx2) const {
        unsigned size = generators.size();
        auto permuted = generators;

        for (auto &gen : permuted)
            std::swap(gen.eq.a[idx1], gen.eq.a[idx2]);

        std::vector<bool> has_partner(size, false);
        unsigned partner_count = 0;

        for (const auto &gen : generators) {
            unsigned i;
            for (i = 0; i < size; ++i) {
                if (has_partner[i])
                    continue;

                if (gen.eq == permuted[i].eq) {
                    has_partner[i] = true;
                    ++partner_count;
                    break;
                }
            }
            if (i == size)
                return false;
        }

        return partner_count == size;
    }

    // the indices in the pairs are in ascending orders
    __attribute__ ((noinline))
    std::vector<std::pair<unsigned, unsigned>> symmetries() const {
        std::vector<std::pair<unsigned, unsigned>> ret{};

        std::array<bool, dim> done{};
        done.fill(false);

        for (unsigned i = 0; i < dim - 1; ++i) {
            if (done[i])
                continue;

            unsigned k = i;
            for (unsigned j = i + 1; j < dim; ++j) {
                if (done[j])
                    continue;

                if (has_symmetry(i, j)) {
                    ret.push_back(std::pair<unsigned, unsigned>(k, j));
                    k = j;
                    done[j] = true;
                }
            }
        }

        return ret;
    }

    __attribute__ ((noinline))
    WeightSystemBuilder restrict(const Vector &x) const {
        WeightSystemBuilder ret{};
        ret.iteration_ = iteration_ + 1;

        // for (int i = 0; i < num_symmetries; ++i) {
        //     if (x[symmetries[i].first] == x[symmetries[i].second]) {
        //         // still a symmetry
        //         ret.symmetries[ret.num_symmetries] = symmetries[i];
        //         ++ret.num_symmetries;
        //         continue;
        //     }
        // }

        std::vector<Long> xq;
        xq.resize(generators.size());

        for (unsigned i = 0; i < generators.size(); i++) {
            xq[i] = generators[i].eq.apply_to(x);
            if (xq[i] == 0)
                ret.generators.push_back(generators[i]);
        }

        for (unsigned i = 0; i < generators.size() - 1; i++) {
            auto &c1 = generators[i];

            for (unsigned j = i + 1; j < generators.size(); j++) {
                auto &c2 = generators[j];

                if (xq[i] * xq[j] >= 0)
                    continue;

                auto new_incidences = c1.incidences | c2.incidences;
                if (new_incidences.count() > iteration_ + 2)
                    continue;

                bool is_redundant =
                    std::any_of(ret.generators.begin(), ret.generators.end(),
                                [&](Generator &c) {
                                    return c.incidences < new_incidences;
                                });
                if (is_redundant)
                    continue;

                rearranging_erase_if(ret.generators, [&](Generator &c) {
                    return new_incidences < c.incidences;
                });

                auto new_eq = intersect(c1.eq, c2.eq, x);
                if (new_eq.c < 0)
                    new_eq = -new_eq;

                ret.generators.push_back(Generator{new_eq, new_incidences});
            }
        }

        return ret;
    }

    __attribute__ ((noinline))
    bool average_if_nonzero(WeightSystem &q) const {
        q.c = 1;
        for (size_t i = 0; i < generators.size(); ++i)
            q.c = std::experimental::lcm(q.c, generators[i].eq.c);

        for (size_t j = 0; j < dim; ++j) {
            q.a[j] = 0;
            for (size_t i = 0; i < generators.size(); ++i)
                q.a[j] += generators[i].eq.a[j] * (q.c / generators[i].eq.c);

            if (q.a[j] == 0)
                return false;
        }

        q.c *= static_cast<Long>(generators.size());
        q.cancel();

        return true;
    }

    friend std::ostream &operator<<(std::ostream &os,
                                    const WeightSystemBuilder &rhs) {
        for (auto &generator : rhs.generators)
            os << generator.eq << std::endl;
        return os;
    }

    __attribute__ ((noinline))
    static bool gives_good_weightsystem(const Vector &x, Long r_numerator,
                                        Long r_denominator) {
        Long xsum = 0, xmax = 0;

        for (int l = 0; l < dim; l++) {
            xsum += x[l];
            if (x[l] > xmax)
                xmax = x[l];
        }

        // Point leads to weight systems containing a weight of 1
        if (xsum < 2)
            return false;

        // Point leads to weight systems containing a weight of 1/2 or two
        // weights with a sum of 1
        if (xsum == 2)
            return false;

        // Point does not allow positive weight systems (except if all
        // coordinates are 1/r)
        if (xmax * r_numerator <= r_denominator)
            return false;

        return true;
    }

private:
    WeightSystemBuilder() {}
};

#endif

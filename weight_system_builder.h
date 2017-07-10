#ifndef WEIGHT_SYSTEM_BUILDER
#define WEIGHT_SYSTEM_BUILDER

#include <algorithm>
#include <bitset>
#include <experimental/numeric>
#include <vector>
#include "stl_utils.h"
#include "vector.h"
#include "weight_system.h"

class WeightSystemBuilder {
    struct Generator {
        WeightSystem eq;
        std::bitset<dim> incidences;
    };
    std::vector<Generator>
        generators; // TODO: can there ever be more generators than dimensions?
    unsigned iteration_;

public:
    // bool operator<(const WeightSystemBuilder &rhs) const
    // {
    //     if (generators.size() != rhs.generators.size())
    //         return generators.size() < rhs.generators.size();

    //     for (unsigned i = 0; i < generators.size(); ++i)
    //         if (generators[i].eq != rhs.generators[i].eq)
    //             return generators[i].eq < rhs.generators[i].eq;

    //     return false;
    // }

    // private

    // WeightSystemBuilder canonicalized() {
    //     using std::array;
    //     using std::vector;

    //     unsigned size = generators.size();

    //     vector<array<vector<Long>, dim>> gens_perms;

    //     for (unsigned i = 0; i < size; ++i) {
    //         array<vector<Long>, dim> gens{};
    //         gens.

    //     }

    //     for (unsigned i = 0; i < size; ++i)
    //         std::sort(Muh{}, Muh{}, [](auto &a, auto &b) { return true; });

    //     return *this;
    // }

    WeightSystemBuilder() : iteration_{0}
    {
        // num_symmetries = dim - 1;
        // for (int i = 0; i < dim - 1; ++i)
        //     symmetries[i] = {i, i + 1};

        generators.resize(dim);
        for (int generator_nr = 0; generator_nr < dim; ++generator_nr) {
            auto &generator = generators[generator_nr];

            for (int i = 0; i < dim; ++i)
                generator.eq.weights[i] = generator_nr == i ? 1 : 0;

            generator.incidences.set(generator_nr);
        }
    }

    unsigned iteration() const { return iteration_; }
    unsigned generator_count() const
    {
        return static_cast<unsigned>(generators.size());
    }

    bool has_symmetry(unsigned idx1, unsigned idx2) const
    {
        size_t size = generators.size();
        auto permuted = generators;

        for (auto &gen : permuted)
            std::swap(gen.eq.weights[idx1], gen.eq.weights[idx2]);

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
    __attribute__((noinline)) std::vector<std::pair<unsigned, unsigned>>
    symmetries() const
    {
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

    __attribute__((noinline)) WeightSystemBuilder restrict(
        const Vector &x) const
    {
        WeightSystemBuilder ret{Noinit{}};
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
            xq[i] = distance(generators[i].eq, x);
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
                // if (new_eq.c < 0)
                //     new_eq = -new_eq;

                ret.generators.push_back(Generator{new_eq, new_incidences});
            }
        }

        return ret;
    }

    // __attribute__ ((noinline))
    // bool sum_if_nonzero(WeightSystem &q) const {
    //     for (unsigned j = 0; j < dim; ++j) {
    //         q.weights[j] = 0;
    //         for (unsigned i = 0; i < generators.size(); ++i)
    //             q.weights[j] += generators[i].eq.weights[j];

    //         if (q.weights[j] == 0)
    //             return false;
    //     }

    //     q.cancel();

    //     return true;
    // }

    __attribute__((noinline)) bool average_if_nonzero(WeightSystem &q) const
    {
        size_t size = generators.size();
        if (size == 0)
            return false;

        std::vector<Long> norms;
        norms.reserve(size);

        for (unsigned i = 0; i < size; ++i)
            norms.push_back(norm(generators[i].eq));

        Long lcm = norms[0];
        for (unsigned i = 1; i < size; ++i)
            lcm = std::experimental::lcm(lcm, norms[i]);

        for (unsigned j = 0; j < dim; ++j) {
            q.weights[j] = 0;
            for (unsigned i = 0; i < size; ++i)
                q.weights[j] += generators[i].eq.weights[j] * (lcm / norms[i]);

            if (q.weights[j] == 0)
                return false;
        }

        cancel(q);

        return true;
    }

    friend std::ostream &operator<<(std::ostream &os,
                                    const WeightSystemBuilder &rhs)
    {
        for (auto &generator : rhs.generators)
            os << generator.eq << std::endl;
        return os;
    }

    bool allows(const Vector &x) const
    {
        for (const auto &gen : generators)
            if (distance(gen.eq, x) != 0)
                return false;
        return true;
    }

    __attribute__((noinline)) static bool leads_to_allowed_weightsystem(
        const Vector &x, Long r_numerator, Long r_denominator)
    {
        Long xsum = 0, xmax = 0;

        for (int l = 0; l < dim; l++) {
            xsum += x.coords[l];
            if (x.coords[l] > xmax)
                xmax = x.coords[l];
        }

        // Point is not allowed by a positive weight system
        if (xmax * r_numerator <= r_denominator)
            return false;

        if (!allow_weight_one && xsum == 1)
            return false;

        if (!allow_weight_one_half && xsum == 2 && xmax == 2)
            return false;

        if (!allow_weights_sum_one && xsum == 2 && xmax != 2)
            return false;

        return true;
    }

private:
    struct Noinit {
    };

    WeightSystemBuilder(Noinit) {}
};

#endif

#include "weight_system_builder.h"
#include <algorithm>
#include <experimental/numeric>
#include "stl_utils.h"

WeightSystemBuilder::WeightSystemBuilder() : iteration{0}
{
    generators.resize(dim);
    for (int generator_nr = 0; generator_nr < dim; ++generator_nr) {
        auto &generator = generators[generator_nr];

        for (int i = 0; i < dim; ++i)
            generator.eq.weights[i] = generator_nr == i ? 1 : 0;

        generator.incidences.set(generator_nr);
    }
}

bool WeightSystemBuilder::operator<(const WeightSystemBuilder &rhs) const
{
    if (generators.size() != rhs.generators.size())
        return generators.size() < rhs.generators.size();

    for (unsigned i = 0; i < generators.size(); ++i)
        if (generators[i].eq != rhs.generators[i].eq)
            return generators[i].eq < rhs.generators[i].eq;

    return false;
}

unsigned WeightSystemBuilder::generator_count() const
{
    return static_cast<unsigned>(generators.size());
}

bool WeightSystemBuilder::has_symmetry(unsigned idx1, unsigned idx2) const
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
std::vector<std::pair<unsigned, unsigned>> WeightSystemBuilder::symmetries()
    const
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

WeightSystemBuilder WeightSystemBuilder::restrict(const Vector &x) const
{
    WeightSystemBuilder ret{NoInit{}};
    ret.iteration = iteration + 1;

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
            if (new_incidences.count() > iteration + 2)
                continue;

            bool is_redundant = std::any_of(
                ret.generators.begin(), ret.generators.end(),
                [&](Generator &c) { return c.incidences < new_incidences; });
            if (is_redundant)
                continue;

            rearranging_erase_if(ret.generators, [&](Generator &c) {
                return new_incidences < c.incidences;
            });

            auto new_eq = intersect(c1.eq, c2.eq, x);
            ret.generators.push_back(Generator{new_eq, new_incidences});
        }
    }

    return ret;
}

bool WeightSystemBuilder::average_if_nonzero(WeightSystem &q) const
{
    size_t size = generators.size();
    if (size == 0)
        return false;

    std::vector<Long> norms;
    norms.reserve(size);

    for (const auto &gen : generators)
        norms.push_back(norm(gen.eq));

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

std::ostream &operator<<(std::ostream &os, const WeightSystemBuilder &rhs)
{
    for (auto &generator : rhs.generators)
        os << generator.eq << std::endl;
    return os;
}

bool WeightSystemBuilder::allows(const Vector &x) const
{
    for (const auto &gen : generators)
        if (distance(gen.eq, x) != 0)
            return false;
    return true;
}

bool leads_to_allowed_weightsystem(const Vector &x)
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

#include "polytope.h"
#include "palp/Global.h"

using std::endl;

int euler_number(const PolytopeInfo &info)
{
    assert(dim == 6);
    return 48 +
           6 * (info.hodge_numbers_1[1] - info.hodge_numbers_1[2] +
                info.hodge_numbers_1[3]);
}

std::ostream &operator<<(std::ostream &os, const PolytopeInfo &info)
{
    if (!info.ip)
        return os << "not ip";

    os << "M:" << info.point_count << " " << info.vertex_count;

    if (info.reflexive) {
        os << " N:" << info.dual_point_count << " " << info.facet_count;
        os << " H:" << info.hodge_numbers_1[1];
        for (unsigned i = 2; i < info.hodge_numbers_1.size(); ++i)
            os << "," << info.hodge_numbers_1[i];

        if (dim == 6)
            os << " [" << euler_number(info) << "]";
    } else {
        os << " F:" << info.facet_count;
    }

    return os;
}

void read(BufferedReader &f, PolytopeInfo &info)
{
    info = PolytopeInfo{};

    uint8_t type{};
    read(f, type);

    switch (type) {
    case 0:
        info.ip = false;
        return;
    case 1:
        info.ip = true;
        info.reflexive = false;
        break;
    case 2:
        info.ip = true;
        info.reflexive = true;
        break;
    }

    info.vertex_count = static_cast<unsigned>(read_varint(f));
    info.facet_count = static_cast<unsigned>(read_varint(f));
    info.point_count = static_cast<unsigned>(read_varint(f));

    if (info.reflexive) {
        info.dual_point_count = static_cast<unsigned>(read_varint(f));
        for (unsigned i = 1; i < info.hodge_numbers_1.size(); ++i)
            info.hodge_numbers_1[i] = static_cast<unsigned>(read_varint(f));
    }
}

void write(BufferedWriter &f, const PolytopeInfo &info)
{
    if (!info.ip) {
        write(f, static_cast<uint8_t>(0));
        return;
    }

    write(f, static_cast<uint8_t>(info.reflexive ? 2 : 1));

    write_varint(f, info.vertex_count);
    write_varint(f, info.facet_count);
    write_varint(f, info.point_count);

    if (info.reflexive) {
        write_varint(f, info.dual_point_count);
        for (unsigned i = 1; i < info.hodge_numbers_1.size(); ++i)
            write_varint(f, info.hodge_numbers_1[i]);
    }
}

extern "C" int QuickAnalysis(PolyPointList *_P, BaHo *_BH, FaceInfo *_FI);

template <typename T>
static int log2(T w)
{
    int i = -1;
    while (w) {
        w /= 2;
        i++;
    }
    return i;
}

template <typename T>
static void update_max(T &max, const T &current)
{
    if (current > max)
        max = current;
}

template <typename T>
static void update_min(T &min, const T &current)
{
    if (current < min)
        min = current;
}

void analyze(const WeightSystem<dim> &ws, PolytopeInfo &info,
             PolytopeStatistics &stats)
{
    using std::make_unique;

    assert(r_denominator == 1);
    assert(POLY_Dmax >= dim - 1);

    // restrictions of the current implementation
    assert(dim == 6);
    assert(r_numerator == 1);

    // allocate large structures on the heap once
    static auto CW = make_unique<CWS>();
    static auto P = make_unique<PolyPointList>();
    static auto FI = make_unique<FaceInfo>();

    Ring greatest_weight = 0;

    CW->N = dim;
    CW->index = r_numerator;
    CW->nw = 1;
    Ring norm = 0;
    for (unsigned i = 0; i < dim; i++) {
        CW->W[0][i] = ws.weights[i] * r_numerator;
        norm += ws.weights[i];
        update_max(greatest_weight, ws.weights[i]);
    }
    CW->d[0] = norm;
    CW->nz = 0;

    Make_CWS_Points(&*CW, &*P);

    BaHo BH;

    info = PolytopeInfo{};

    info.ip = QuickAnalysis(&*P, &BH, &*FI);

    if (!info.ip) {
        stats.n_nonIP++;
        return;
    }

    assert(BH.n == dim - 1);

    info.reflexive = BH.np > 0;

    int log_w = log2(CW->W[0][5]);
    assert(0 <= log_w && log_w < PolytopeStatistics::max_log2);

    info.vertex_count = BH.mv;
    info.facet_count = BH.nv;
    info.point_count = BH.mp;

    if (info.reflexive) {
        stats.n_ref++;
        stats.n_w[log_w]++;

        info.dual_point_count = BH.np;

        for (unsigned i = 1; i < info.hodge_numbers_1.size(); ++i)
            info.hodge_numbers_1[i] = BH.h1[i];

        int chi = euler_number(info);

        update_max(stats.max_mp, BH.mp);
        update_max(stats.max_mv, BH.mv);
        update_max(stats.max_np, BH.np);
        update_max(stats.max_nv, BH.nv);
        update_max(stats.max_h22, BH.h22);
        update_max(stats.max_chi, chi);
        update_min(stats.min_chi, chi);
        update_max(stats.max_w, greatest_weight);

        for (unsigned i = 1; i < stats.max_h1.size(); ++i)
            update_max(stats.max_h1[i], BH.h1[i]);
        for (unsigned i = 0; i < stats.max_nf.size(); ++i)
            update_max(stats.max_nf[i], FI->nf[i]);
    } else {
        stats.n_IP_nonRef++;
        stats.nr_n_w[log_w]++;

        update_max(stats.nr_max_mp, BH.mp);
        update_max(stats.nr_max_mv, BH.mv);
        update_max(stats.nr_max_nv, BH.nv);
        update_max(stats.nr_max_w, greatest_weight);
    }
}

std::ostream &operator<<(std::ostream &os, const PolytopeStatistics &stats)
{
    os << "non-IP: #=" << stats.n_nonIP << endl;
    os << "IP, non-reflexive: #=" << stats.n_IP_nonRef
       << ", max_mp=" << stats.nr_max_mp << ", max_mv=" << stats.nr_max_mv
       << ", max_nv=" << stats.nr_max_nv << ", max_w=" << stats.nr_max_w
       << endl;

    os << "  #(greatest weight) of given log2: ";
    for (unsigned i = 0; i < stats.nr_n_w.size(); i++)
        os << " " << i << ":" << stats.nr_n_w[i];
    os << endl;

    os << "reflexive: #=" << stats.n_ref << ", max_mp=" << stats.max_mp
       << ", max_mv=" << stats.max_mv << ", max_np=" << stats.max_np
       << ", max_nv=" << stats.max_nv << ", max_w=" << stats.max_w << endl;

    os << "  #(greatest weight) of given log2: ";
    for (unsigned i = 0; i < stats.n_w.size(); ++i)
        os << " " << i << ":" << stats.n_w[i];
    os << endl;

    os << "  max #(faces): ";
    for (auto nf : stats.max_nf)
        os << nf << " ";
    os << endl;

    os << "  ";
    for (unsigned i = 1; i < stats.max_h1.size(); ++i)
        os << "h1" << i << "<=" << stats.max_h1[i] << ", ";
    os << "h22<=" << stats.max_h22 << ", ";

    os << stats.min_chi << "<=chi<=" << stats.max_chi << endl;

    return os;
}

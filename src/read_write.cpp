#include "read_write.h"
#include <cassert>

void write_config(BufferedWriter &f, unsigned dim1, Ring r_numerator1,
                  Ring r_denominator1)
{
    write32u(f, dim1);
    write32u(f, r_numerator1);
    write32u(f, r_denominator1);
}

void write_config(BufferedWriter &f)
{
    write_config(f, dim, r_numerator, r_denominator);
}

void check_config(BufferedReader &f, unsigned dim1, Ring r_numerator1,
                  Ring r_denominator1)
{
    uint32_t u;
    int32_t i;

    read(f, u);
    assert(u == dim1);
    read(f, i);
    assert(i == r_numerator1);
    read(f, i);
    assert(i == r_denominator1);
}

void check_config(BufferedReader &f)
{
    check_config(f, dim, r_numerator, r_denominator);
}

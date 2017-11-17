#ifndef READ_WRITE_H
#define READ_WRITE_H

#include "buffered_reader.h"
#include "buffered_writer.h"
#include "config.h"

void write_config(BufferedWriter &f);
void write_config(BufferedWriter &f, unsigned dim1, Ring r_numerator1,
                  Ring r_denominator1);

void check_config(BufferedReader &f);
void check_config(BufferedReader &f, unsigned dim1, Ring r_numerator1,
                  Ring r_denominator1);

#endif

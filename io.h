#ifndef IO_H
#define IO_H

#include <iostream>

void read(std::istream &f, uint16_t &data);
void read(std::istream &f, uint32_t &data);
void read(std::istream &f, int16_t &data);
void read(std::istream &f, int32_t &data);
void write(std::ostream &f, uint16_t data);
void write(std::ostream &f, uint32_t data);
void write(std::ostream &f, int16_t data);
void write(std::ostream &f, int32_t data);

#endif

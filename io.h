#ifndef IO_H
#define IO_H

#include <boost/optional.hpp>
#include <exception>
#include <iostream>

struct VarintReadError : std::exception {
    virtual const char *what() const noexcept { return "Varint read error"; }
};

void read(std::istream &f, uint16_t &data);
void read(std::istream &f, uint32_t &data);
void read(std::istream &f, int16_t &data);
void read(std::istream &f, int32_t &data);
void write(std::ostream &f, uint16_t data);
void write(std::ostream &f, uint32_t data);
void write(std::ostream &f, int16_t data);
void write(std::ostream &f, int32_t data);

void write_varint(std::ostream &f, unsigned long i);
unsigned long read_varint(std::istream &f);

#endif

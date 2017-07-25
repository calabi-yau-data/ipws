#include "io.h"
#include <arpa/inet.h>

void read(std::istream &f, uint16_t &data)
{
    f.read(reinterpret_cast<char *>(&data), sizeof(data));
    data = ntohs(data);
}

void read(std::istream &f, uint32_t &data)
{
    f.read(reinterpret_cast<char *>(&data), sizeof(data));
    data = ntohl(data);
}

void read(std::istream &f, int16_t &data)
{
    read(f, *reinterpret_cast<uint16_t *>(&data));
}

void read(std::istream &f, int32_t &data)
{
    read(f, *reinterpret_cast<uint32_t *>(&data));
}

void write(std::ostream &f, uint16_t data)
{
    data = htons(data);
    f.write(reinterpret_cast<const char *>(&data), sizeof(data));
}

void write(std::ostream &f, uint32_t data)
{
    data = htonl(data);
    f.write(reinterpret_cast<const char *>(&data), sizeof(data));
}

void write(std::ostream &f, int16_t data)
{
    write(f, static_cast<uint16_t>(data));
}

void write(std::ostream &f, int32_t data)
{
    write(f, static_cast<uint32_t>(data));
}

void write_varint(std::ostream &f, unsigned long i)
{
    while (i > 127) {
        unsigned char v = static_cast<unsigned char>(i) | 128;
        f.write(reinterpret_cast<char *>(&v), sizeof(v));
        i >>= 7;
    }

    unsigned char v = static_cast<unsigned char>(i);
    f.write(reinterpret_cast<char *>(&v), sizeof(v));
}

unsigned long read_varint(std::istream &f)
{
    unsigned long ret = 0;
    unsigned pos = 0;

    while (true) {
        unsigned char v;
        f.read(reinterpret_cast<char *>(&v), sizeof(v));

        unsigned long w = v & 127;

        ret |= w << pos;

        if (ret >> pos != w)
            throw VarintReadError();

        if ((v & 128) == 0)
            break;

        pos += 7;
    }

    return ret;
}

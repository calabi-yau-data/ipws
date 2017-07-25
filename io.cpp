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

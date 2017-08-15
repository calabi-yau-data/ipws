#include "buffered_reader.h"
#include <endian.h>
#include <memory.h>
#include <algorithm>

const char *BufferedReader::EofError::what() const noexcept
{
    return "End of file reached while reading";
}

const char *BufferedReader::OpenError::what() const noexcept
{
    return "Could not open file for reading";
}

BufferedReader::BufferedReader(const std::string &path)
        : stream{}, buffer{}, buffer_data_start{0}, buffer_data_size{0}
{
    buffer.resize(100000);

    stream.open(path, std::ios::binary);
    if (!stream)
        throw OpenError{};
    stream.exceptions(std::ios::badbit);
}

void BufferedReader::read(void *data, size_t size)
{
    char *dest = reinterpret_cast<char *>(data);

    while (true) {
        size_t copy_now = std::min(size, buffer_data_size);

        // std::copy(buffer.data() + buffer_data_start,
        //           buffer.data() + buffer_data_start + copy_now, dest);
        memcpy(dest, buffer.data() + buffer_data_start, copy_now);
        dest += copy_now;
        size -= copy_now;

        if (size == 0) {
            buffer_data_start += copy_now;
            buffer_data_size -= copy_now;
            return;
        }

        stream.read(buffer.data(), buffer.size());

        buffer_data_start = 0;
        buffer_data_size = stream.gcount();

        if (buffer_data_size == 0)
            throw EofError{};
    }
}

void BufferedReader::seek(size_t pos)
{
    buffer_data_size = 0;
    stream.seekg(pos);
}

void BufferedReader::seek_relative(ptrdiff_t offset)
{
    buffer_data_size = 0;
    stream.seekg(offset, std::ios_base::cur);
}

void read(BufferedReader &f, uint16_t &data)
{
    f.read(&data, sizeof(data));
    data = be16toh(data);
}

void read(BufferedReader &f, uint32_t &data)
{
    f.read(&data, sizeof(data));
    data = be32toh(data);
}

void read(BufferedReader &f, uint64_t &data)
{
    f.read(&data, sizeof(data));
    data = be64toh(data);
}

void read(BufferedReader &f, int16_t &data)
{
    read(f, *reinterpret_cast<uint16_t *>(&data));
}

void read(BufferedReader &f, int32_t &data)
{
    read(f, *reinterpret_cast<uint32_t *>(&data));
}

void read(BufferedReader &f, int64_t &data)
{
    read(f, *reinterpret_cast<uint64_t *>(&data));
}

unsigned long read_varint(BufferedReader &f)
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

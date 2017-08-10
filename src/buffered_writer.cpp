#include "buffered_writer.h"
#include <arpa/inet.h>
#include <memory.h>
#include <algorithm>

const char *BufferedWriter::OpenError::what() const noexcept
{
    return "Could not open file for writing";
}

BufferedWriter::BufferedWriter(const std::string &path)
        : stream{}, buffer{}, buffer_data_end{0}
{
    stream.open(path, std::ios::binary);

    if (!stream)
        throw OpenError{};

    stream.exceptions(std::ios::badbit);
}

BufferedWriter::~BufferedWriter()
{
    flush();
}

void BufferedWriter::write(const void *data, unsigned size)
{
    const char *src = reinterpret_cast<const char *>(data);

    while (true) {
        unsigned copy_now = std::min(
            size, static_cast<unsigned>(buffer.size()) - buffer_data_end);

        // std::copy(src, src + copy_now, buffer.data() + buffer_data_end);
        memcpy(buffer.data() + buffer_data_end, src, copy_now);

        src += copy_now;
        size -= copy_now;
        buffer_data_end += copy_now;

        if (buffer_data_end != buffer.size())
            return;

        flush();
    }
}

void BufferedWriter::flush()
{
    flush_to_stream();
    stream.flush();
}

void BufferedWriter::flush_to_stream()
{
    stream.write(buffer.data(), buffer_data_end);
    buffer_data_end = 0;
}

void BufferedWriter::seek(unsigned pos)
{
    flush_to_stream();
    stream.seekp(pos);
}

void BufferedWriter::seek_relative(int offset)
{
    flush_to_stream();
    stream.seekp(offset, std::ios_base::cur);
}

void write(BufferedWriter &f, uint16_t data)
{
    data = htons(data);
    f.write(&data, sizeof(data));
}

void write(BufferedWriter &f, uint32_t data)
{
    data = htonl(data);
    f.write(&data, sizeof(data));
}

void write(BufferedWriter &f, int16_t data)
{
    write(f, static_cast<uint16_t>(data));
}

void write(BufferedWriter &f, int32_t data)
{
    write(f, static_cast<uint32_t>(data));
}

void write_varint(BufferedWriter &f, unsigned long i)
{
    while (i > 127) {
        unsigned char v = static_cast<unsigned char>(i) | 128;
        f.write(reinterpret_cast<char *>(&v), sizeof(v));
        i >>= 7;
    }

    unsigned char v = static_cast<unsigned char>(i);
    f.write(reinterpret_cast<char *>(&v), sizeof(v));
}
#include "buffered_writer.h"
#include <algorithm>
#include <endian.h>
#include <memory.h>

const char *BufferedWriter::OpenError::what() const noexcept
{
    return "Could not open file for writing";
}

BufferedWriter::BufferedWriter(const std::string &path, size_t buffer_size)
        : stream{}, buffer{}, buffer_data_end{0}
{
    buffer.resize(buffer_size);

    stream.open(path, std::ios::binary);
    if (!stream)
        throw OpenError{};
    stream.exceptions(std::ios::badbit);
}

BufferedWriter::~BufferedWriter()
{
    flush();
}

void BufferedWriter::write(const void *data, size_t size)
{
    const char *src = reinterpret_cast<const char *>(data);

    while (true) {
        size_t copy_now = std::min(size, buffer.size() - buffer_data_end);

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

void BufferedWriter::seek(size_t pos)
{
    flush_to_stream();
    stream.seekp(pos);
}

void BufferedWriter::seek_relative(std::ptrdiff_t offset)
{
    flush_to_stream();
    stream.seekp(offset, std::ios_base::cur);
}

void write8u(BufferedWriter &f, uint8_t data)
{
    f.write(&data, sizeof(data));
}

void write16u(BufferedWriter &f, uint16_t data)
{
    data = htobe16(data);
    f.write(&data, sizeof(data));
}

void write32u(BufferedWriter &f, uint32_t data)
{
    data = htobe32(data);
    f.write(&data, sizeof(data));
}

void write64u(BufferedWriter &f, uint64_t data)
{
    data = htobe64(data);
    f.write(&data, sizeof(data));
}

void write8i(BufferedWriter &f, int8_t data)
{
    write8u(f, static_cast<uint8_t>(data));
}

void write16i(BufferedWriter &f, int16_t data)
{
    write16u(f, static_cast<uint16_t>(data));
}

void write32i(BufferedWriter &f, int32_t data)
{
    write32u(f, static_cast<uint32_t>(data));
}

void write64i(BufferedWriter &f, int64_t data)
{
    write64u(f, static_cast<uint64_t>(data));
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

unsigned varint_storage_size(unsigned long i)
{
    unsigned size = 1;
    while (i > 127) {
        ++size;
        i >>= 7;
    }
    return size;
}

#ifndef BUFFERED_WRITER_H
#define BUFFERED_WRITER_H

#include <exception>
#include <fstream>
#include <vector>

class BufferedWriter {
public:
    BufferedWriter(const std::string &path, size_t buffer_size = 100000);
    ~BufferedWriter();

    BufferedWriter(const BufferedWriter &other) = delete;
    BufferedWriter &operator=(const BufferedWriter &rhs) = delete;

    BufferedWriter(BufferedWriter &&other) = default;
    BufferedWriter &operator=(BufferedWriter &&rhs) = default;

    void seek(size_t pos);
    void seek_relative(std::ptrdiff_t offset);

    void write(const void *data, size_t size);
    void flush();

    struct OpenError : std::exception {
        virtual const char *what() const noexcept;
    };

private:
    void flush_to_stream();

    std::ofstream stream;
    std::vector<char> buffer;
    size_t buffer_data_end;
};

void write(BufferedWriter &f, uint8_t data);
void write(BufferedWriter &f, uint16_t data);
void write(BufferedWriter &f, uint32_t data);
void write(BufferedWriter &f, uint64_t data);
void write(BufferedWriter &f, int8_t data);
void write(BufferedWriter &f, int16_t data);
void write(BufferedWriter &f, int32_t data);
void write(BufferedWriter &f, int64_t data);

void write_varint(BufferedWriter &f, unsigned long i);

#endif

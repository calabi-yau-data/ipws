#ifndef BUFFERED_WRITER_H
#define BUFFERED_WRITER_H

#include <array>
#include <exception>
#include <fstream>

class BufferedWriter {
public:
    BufferedWriter(const std::string &path);
    ~BufferedWriter();

    BufferedWriter(const BufferedWriter &other) = delete;
    BufferedWriter &operator=(const BufferedWriter &rhs) = delete;

    BufferedWriter(BufferedWriter &&other) = default;
    BufferedWriter &operator=(BufferedWriter &&rhs) = default;

    void seek(unsigned pos);
    void seek_relative(int offset);

    void write(const void *data, unsigned size);
    void flush();

    struct OpenError : std::exception {
        virtual const char *what() const noexcept;
    };

private:
    void flush_to_stream();

    std::ofstream stream;
    std::array<char, 1024> buffer;
    unsigned buffer_data_end;
};

void write(BufferedWriter &f, uint16_t data);
void write(BufferedWriter &f, uint32_t data);
void write(BufferedWriter &f, int16_t data);
void write(BufferedWriter &f, int32_t data);

void write_varint(BufferedWriter &f, unsigned long i);

#endif

#ifndef BUFFERED_READER_H
#define BUFFERED_READER_H

#include <exception>
#include <fstream>
#include <vector>

struct VarintReadError : std::exception {
    virtual const char *what() const noexcept { return "Varint read error"; }
};

class BufferedReader {
public:
    BufferedReader(const std::string &path);

    BufferedReader(const BufferedReader &other) = delete;
    BufferedReader &operator=(const BufferedReader &rhs) = delete;

    BufferedReader(BufferedReader &&other) = default;
    BufferedReader &operator=(BufferedReader &&rhs) = default;

    void seek(size_t pos);
    void seek_relative(ptrdiff_t offset);

    void read(void *data, size_t size);

    struct EofError : std::exception {
        virtual const char *what() const noexcept;
    };

    struct OpenError : std::exception {
        virtual const char *what() const noexcept;
    };

private:
    std::ifstream stream;
    std::vector<char> buffer;
    size_t buffer_data_start;
    size_t buffer_data_size;
};

void read(BufferedReader &f, uint8_t &data);
void read(BufferedReader &f, uint16_t &data);
void read(BufferedReader &f, uint32_t &data);
void read(BufferedReader &f, uint64_t &data);
void read(BufferedReader &f, int8_t &data);
void read(BufferedReader &f, int16_t &data);
void read(BufferedReader &f, int32_t &data);
void read(BufferedReader &f, int64_t &data);

unsigned long read_varint(BufferedReader &f);

#endif

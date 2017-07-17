#ifndef FILE_H
#define FILE_H

#include <experimental/optional>
#include <gsl/gsl>

class File {
public:
    template <class T>
    using optional = std::experimental::optional<T>;

    class Error {
    };

    File();
    File(const File &other);
    File(File &&other);
    ~File();
    File &operator=(File rhs);
    friend void swap(File &a, File &b);

    static optional<File> open(const std::string &path);
    static optional<File> create_new(const std::string &path);

    void seek(unsigned pos);

    void read(void *data, size_t size);
    void read(uint16_t &data);
    void read(uint32_t &data);

    void write(const void *data, size_t size);
    void write(uint16_t data);
    void write(uint32_t data);

    template <class T>
    void read(gsl::span<T> data)
    {
        for (auto &x : data)
            read(x);
    }

    template <class T>
    void write(gsl::span<T> data)
    {
        for (auto &x : data)
            write(x);
    }

private:
    int fd;

    File(int fd);
};

#endif

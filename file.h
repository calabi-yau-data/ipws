#ifndef FILE_H
#define FILE_H

#include <boost/filesystem.hpp>
#include <boost/optional.hpp>
#include <exception>
#include <gsl/gsl>

class File {
public:
    template <class T>
    using optional = boost::optional<T>;

    struct Error : std::exception {
        virtual const char *what() const noexcept { return "IO error"; }
    };

    File();
    File(const File &other);
    File(File &&other);
    ~File();
    File &operator=(File rhs);
    friend void swap(File &a, File &b);

    static optional<File> open(const boost::filesystem::path &path);
    static optional<File> create_new(const boost::filesystem::path &path);

    void seek(unsigned pos);
    void seek_relative(int pos);

    void read(void *data, size_t size);
    void read(uint16_t &data);
    void read(uint32_t &data);
    void read(int16_t &data);
    void read(int32_t &data);

    void write(const void *data, size_t size);
    void write(uint16_t data);
    void write(uint32_t data);
    void write(int16_t data);
    void write(int32_t data);

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

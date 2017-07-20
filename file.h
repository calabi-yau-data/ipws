#ifndef FILE_H
#define FILE_H

#include <boost/filesystem.hpp>
#include <boost/optional.hpp>
#include <exception>
#include <gsl/gsl>
#include <memory>

class File {
public:
    template <class T>
    using optional = boost::optional<T>;

    struct Error : std::exception {
        virtual const char *what() const noexcept { return "IO error"; }
    };

    File() = delete;

    File(const File &other) = default;
    File(File &&other) = default;
    File &operator=(const File &rhs) = default;
    File &operator=(File &&rhs) = default;

    static optional<File> open(const boost::filesystem::path &path);
    static optional<File> create_new(const boost::filesystem::path &path);

    void seek(unsigned pos);
    void seek_relative(int pos);

    void read(void *data, size_t size);
    void write(const void *data, size_t size);

private:
    struct Impl;
    std::shared_ptr<Impl> impl;

    File(int fd);
};

void read(File &f, uint16_t &data);
void read(File &f, uint32_t &data);
void read(File &f, int16_t &data);
void read(File &f, int32_t &data);
void write(File &f, uint16_t data);
void write(File &f, uint32_t data);
void write(File &f, int16_t data);
void write(File &f, int32_t data);

#endif

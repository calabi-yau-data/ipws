#include "file.h"
#include <arpa/inet.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

using boost::optional;

File::File() : fd{-1}
{
}

File::File(const File &other)
{
    fd = dup(other.fd);
}

File::File(File &&other) : fd{-1}
{
    swap(*this, other);
}

File::~File()
{
    close(fd);
}

File &File::operator=(File rhs)
{
    swap(*this, rhs);
    return *this;
}

optional<File> File::open(const std::string &path)
{
    int fd = ::open(path.data(), O_RDONLY);
    if (fd == -1)
        return {};
    return File{fd};
}

optional<File> File::create_new(const std::string &path)
{
    int fd =
        ::open(path.data(), O_CREAT | O_EXCL | O_NOFOLLOW | O_WRONLY, 0666);
    if (fd == -1)
        return {};
    return File{fd};
}

void File::seek(unsigned pos)
{
    lseek(fd, pos, SEEK_SET);
}

void File::read(void *data, size_t size)
{
    ssize_t ret = ::read(fd, data, size);
    if (ret != static_cast<ssize_t>(size))
        throw Error{};
}

void File::read(uint16_t &data)
{
    read(&data, sizeof data);
    data = ntohs(data);
}

void File::read(uint32_t &data)
{
    read(&data, sizeof data);
    data = ntohl(data);
}

void File::read(int16_t &data)
{
    read(*reinterpret_cast<uint16_t *>(&data));
}

void File::read(int32_t &data)
{
    read(*reinterpret_cast<uint32_t *>(&data));
}

void File::write(const void *data, size_t size)
{
    ssize_t ret = ::write(fd, data, size);
    if (ret != static_cast<ssize_t>(size))
        throw Error{};
}

void File::write(uint16_t data)
{
    data = htons(data);
    write(&data, sizeof data);
}

void File::write(uint32_t data)
{
    data = htonl(data);
    write(&data, sizeof data);
}

void File::write(int16_t data)
{
    write(static_cast<uint16_t>(data));
}

void File::write(int32_t data)
{
    write(static_cast<uint32_t>(data));
}

File::File(int fd) : fd{fd}
{
    this->fd = fd;
}

void swap(File &a, File &b)
{
    std::swap(a.fd, b.fd);
}

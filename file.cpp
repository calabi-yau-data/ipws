#include "file.h"
#include <arpa/inet.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>

namespace ios = boost::iostreams;

using boost::optional;

struct File::Impl {
    ios::file_descriptor ios_fd;
    ios::stream_buffer<ios::file_descriptor> buffer;
    std::iostream stream;

    Impl(int fd)
            : ios_fd(fd, ios::close_handle), buffer(ios_fd), stream(&buffer)
    {
    }
};

optional<File> File::open(const boost::filesystem::path &path)
{
    int fd = ::open(path.c_str(), O_RDONLY);
    if (fd == -1)
        return {};
    return File{fd};
}

optional<File> File::create_new(const boost::filesystem::path &path)
{
    int fd =
        ::open(path.c_str(), O_CREAT | O_EXCL | O_NOFOLLOW | O_WRONLY, 0666);
    if (fd == -1)
        return {};
    return File{fd};
}

void File::seek(unsigned pos)
{
    impl->stream.seekg(pos, std::ios_base::beg);
}

void File::seek_relative(int pos)
{
    impl->stream.seekg(pos, std::ios_base::cur);
}

void File::read(void *data, size_t size)
{
    impl->stream.read(reinterpret_cast<char *>(data), size);
    if (!impl->stream)
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
    impl->stream.write(reinterpret_cast<const char *>(data), size);
    if (!impl->stream)
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

File::File(int fd)
{
    impl = std::make_shared<Impl>(fd);
}

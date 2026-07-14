/*************************************************************************************************
CNFTools -- Copyright (c) 2021, Ashlin Iser, Fahri Taban, KIT - Karlsruhe Institute of Technology

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/

#ifndef SRC_UTIL_STREAMCOMPRESSOR_H_
#define SRC_UTIL_STREAMCOMPRESSOR_H_

#include <archive.h>
#include <archive_entry.h>

#include <array>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <streambuf>
#include <string>

// Supported streaming compression formats (libarchive filters).
enum class CompressionFormat { XZ, GZIP, BZIP2 };

inline const char *compression_suffix(CompressionFormat format)
{
    switch (format)
    {
    case CompressionFormat::GZIP:
        return ".gz";
    case CompressionFormat::BZIP2:
        return ".bz2";
    case CompressionFormat::XZ:
    default:
        return ".xz";
    }
}

class StreamCompressorException : public std::runtime_error
{
public:
    StreamCompressorException() : std::runtime_error("Stream Compressor Exception") {}
    explicit StreamCompressorException(const std::string &msg) : std::runtime_error(msg) {}
    explicit StreamCompressorException(const std::string &msg, archive *arch) : std::runtime_error(msg + ": " + std::string(archive_error_string(arch))) {}
};

class StreamCompressor
{
    unsigned size_;
    unsigned cursor;

    struct archive *arch;
    struct archive_entry *entry;

    int status;
    bool closed;

    void add_filter(CompressionFormat format)
    {
        switch (format)
        {
        case CompressionFormat::GZIP:
            status = archive_write_add_filter_gzip(arch);
            break;
        case CompressionFormat::BZIP2:
            status = archive_write_add_filter_bzip2(arch);
            break;
        case CompressionFormat::XZ:
        default:
            status = archive_write_add_filter_xz(arch);
            break;
        }
        if (status != ARCHIVE_OK)
            throw StreamCompressorException("Error adding compression filter", arch);
    }

public:
    StreamCompressor(const char *output, unsigned size = 0, CompressionFormat format = CompressionFormat::XZ)
        : size_(size), cursor(0), status(0), closed(false)
    {
        arch = archive_write_new();
        status = archive_write_set_format_raw(arch);
        if (status != ARCHIVE_OK)
            throw StreamCompressorException("Error setting format", arch);
        add_filter(format);
        status = archive_write_open_filename(arch, output);
        if (status != ARCHIVE_OK)
            throw StreamCompressorException("Error open archive", arch);

        entry = archive_entry_new();

        std::filesystem::path p(output);
        auto entry_path = p.filename();
        const std::string ext = entry_path.extension().string();
        if (ext == ".xz" || ext == ".gz" || ext == ".bz2")
        {
            entry_path.replace_extension();
        }

        archive_entry_set_pathname(entry, entry_path.c_str());
        if (size != 0)
            resize_entry(size);
        archive_entry_set_filetype(entry, AE_IFREG);
        archive_entry_set_perm(entry, 0644);

        status = archive_write_header(arch, entry);
        if (status != ARCHIVE_OK)
            throw StreamCompressorException("Error writing header", arch);
    }

    ~StreamCompressor()
    {
        if (!closed)
            close();
    }

    void write(const char *buf, unsigned len)
    {
        cursor += len;
        if (size_ != 0 && cursor > size_)
        {
            throw StreamCompressorException("Attempt to write more than announced");
        }

        int bytes_written = archive_write_data(arch, buf, len);
        if (bytes_written != len)
        {
            throw StreamCompressorException("Error writing to archive", arch);
        }
    }

    void resize_entry(size_t size)
    {
        size_ = size;
        archive_entry_set_size(entry, size);
    }

    friend std::istream &operator>>(std::istream &input, StreamCompressor &cmpr)
    {
        input.seekg(0, input.end);
        unsigned length = input.tellg();
        input.seekg(0, input.beg);

        char *buffer = new char[length];
        input.read(buffer, length);
        if (input.fail() && !input.eof())
        {
            throw StreamCompressorException("Error reading from input stream");
        }
        if (archive_entry_size(cmpr.entry) < length)
        {
            cmpr.resize_entry(length);
        }

        cmpr.write(buffer, length);

        delete[] buffer;
        return input;
    }

    void close()
    {
        archive_entry_free(entry);
        status = archive_write_close(arch);
        if (status != ARCHIVE_OK)
        {
            throw StreamCompressorException("Error closing archive", arch);
        }
        status = archive_write_free(arch);
        if (status != ARCHIVE_OK)
        {
            throw StreamCompressorException("Error freeing archive", arch);
        }
        closed = true;
    }
};

// std::streambuf adapter that streams written characters into a StreamCompressor, so that an
// std::ostream (e.g. a redirected std::cout) can write xz-compressed output incrementally
// without buffering the whole payload in memory.
class CompressorStreamBuf : public std::streambuf
{
    StreamCompressor &compressor_;
    std::array<char, 65536> buffer_;

    int flush_buffer()
    {
        const std::ptrdiff_t n = pptr() - pbase();
        if (n > 0)
        {
            compressor_.write(pbase(), static_cast<unsigned>(n));
            pbump(static_cast<int>(-n));
        }
        return 0;
    }

public:
    explicit CompressorStreamBuf(StreamCompressor &compressor) : compressor_(compressor)
    {
        setp(buffer_.data(), buffer_.data() + buffer_.size());
    }

    ~CompressorStreamBuf() override
    {
        try { flush_buffer(); } catch (...) { }
    }

protected:
    int sync() override { return flush_buffer(); }

    int_type overflow(int_type ch) override
    {
        if (flush_buffer() != 0) return traits_type::eof();
        if (!traits_type::eq_int_type(ch, traits_type::eof()))
        {
            *pptr() = traits_type::to_char_type(ch);
            pbump(1);
        }
        return ch;
    }

    std::streamsize xsputn(const char *s, std::streamsize n) override
    {
        if (flush_buffer() != 0) return 0;
        compressor_.write(s, static_cast<unsigned>(n));
        return n;
    }
};

#endif // SRC_UTIL_STREAMCOMPRESSOR_H_

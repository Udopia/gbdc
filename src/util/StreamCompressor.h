/*************************************************************************************************
CNFTools -- Copyright (c) 2021, Markus Iser, KIT - Karlsruhe Institute of Technology

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

#include <iostream>
#include <string>

class StreamCompressorException : public std::runtime_error {
 public:
    StreamCompressorException() : std::runtime_error("Stream Compressor Exception") { }
    explicit StreamCompressorException(const std::string& msg) : std::runtime_error(msg) { }
};

class StreamCompressor {
    unsigned size_;
    unsigned cursor;

    struct archive* a;
    struct archive_entry* entry;

    bool ok;

 public:
    StreamCompressor(const char* output, unsigned size) : size_(size), cursor(0), ok(true) {
        a = archive_write_new();
        // archive_write_add_filter_gzip(a);
        // archive_write_set_format_pax_restricted(a);
        ok &= archive_write_add_filter_xz(a);
        if (!ok) throw StreamCompressorException("Error setting filter");
        ok &= archive_write_open_filename(a, output);
        if (!ok) throw StreamCompressorException("Error open archive");
        entry = archive_entry_new();
        archive_entry_set_pathname(entry, output);
        archive_entry_set_size(entry, size);
        archive_entry_set_filetype(entry, AE_IFREG);
        archive_entry_set_perm(entry, 0644);
        ok &= archive_write_header(a, entry);
        if (!ok) throw StreamCompressorException("Error writing header");
    }

    ~StreamCompressor() {
        archive_entry_free(entry);
        ok &= archive_write_close(a);
        ok &= archive_write_free(a);
        if (!ok) std::cout << "This is not okay!" << std::endl;
    }

    void write(const char* buf, unsigned len) {
        cursor += len;
        if (len > size_) {
            throw StreamCompressorException("Attempt to write more than announced");
        }
        archive_write_data(a, buf, len);
    }

    friend std::istream &operator>>(std::istream& input, StreamCompressor& cmpr) {
        input.seekg(0, input.end);
        unsigned length = input.tellg();
        input.seekg(0, input.beg);

        char* buffer = new char[length];
        input.read(buffer, length);

        cmpr.write(buffer, length);

        delete[] buffer;
        return input;
    }
};

#endif  // SRC_UTIL_STREAMCOMPRESSOR_H_

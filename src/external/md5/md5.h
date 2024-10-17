/** 
 * Markus Iser (KIT): Got this from https://github.com/CommanderBubble/MD5
 * Many Thanks to Michael Lloyd for this neat implemenation of the md5 algorithm.
 * Made some tiny adaptions, wrote a little wrapper (see eof), and copied this from LICENCE:
 * 
 * The MIT License (MIT)
 * Copyright (c) 2015 Michael
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 **/

#ifndef LIB_MD5_MD5_H_
#define LIB_MD5_MD5_H_

/*
 * Size of a standard MD5 signature in bytes.  This definition is for
 * external programs only.  The MD5 routines themselves reference the
 * signature as 4 unsigned 32-bit integers.
 */
const unsigned int MD5_SIZE = (4 * sizeof(unsigned int));   /* 16 */
const unsigned int MD5_STRING_SIZE = 2 * MD5_SIZE + 1;      /* 33 */

namespace md5 {
/*
* The MD5 algorithm works on blocks of characters of 64 bytes.  This
* is an internal value only and is not necessary for external use.
*/
const unsigned int BLOCK_SIZE = 64;

class md5_t {
 public:
    /*
     * md5_t
     *
     * DESCRIPTION:
     * Initialize structure containing state of MD5 computation. (RFC 1321, 3.3: Step 3). 
     * This is for progressive MD5 calculations only. 
     * If you have the complete string available, call it as below.
     * process should be called for each bunch of bytes and after the last
     * process call, finish should be called to get the signature.
     */
    md5_t();

    /*
     * md5_t
     *
     * DESCRIPTION:
     * This function is used to calculate a MD5 signature for a buffer of bytes. 
     * If you only have part of a buffer that you want to process then md5_t, process, and finish should be used.
     *
     * ARGUMENTS:
     * input - A buffer of bytes whose MD5 signature we are calculating.
     * input_length - The length of the buffer.
     * signature_ - A 16 byte buffer that will contain the MD5 signature.
     */
    md5_t(const void* input, const unsigned int input_length, void* signature_ = nullptr);

    /*
     * process
     *
     * DESCRIPTION:
     * This function is used to progressively calculate an MD5 signature some number of bytes at a time.
     *
     * ARGUMENTS:
     * input - A buffer of bytes whose MD5 signature we are calculating.
     * input_length - The length of the buffer.
     */
    void process(const void* input, const unsigned int input_length);

    /*
     * finish
     *
     * DESCRIPTION: 
     * Finish a progressing MD5 calculation and copy the resulting MD5 signature into the result buffer which should be 16 bytes (MD5_SIZE). 
     * After this call, the MD5 structure cannot be used to calculate a new md5, it can only return its signature.
     *
     * ARGUMENTS: 
     * signature_ - A 16 byte buffer that will contain the MD5 signature.
     */
    void finish(void* signature_ = nullptr);

    /*
     * get_sig
     *
     * DESCRIPTION: 
     * Retrieves the previously calculated signature from the MD5 object.
     *
     * ARGUMENTS: 
     * signature_ - A 16 byte buffer that will contain the MD5 signature.
     */
    void get_sig(void* signature_);

    /*
     * get_string
     *
     * DESCRIPTION: Retrieves the previously calculated signature from the MD5 object in printable format.
     *
     * ARGUMENTS: str_ - a string of characters which should be at least 33 bytes long (2 characters per MD5 byte and 1 for the \0).
     */
    void get_string(void* str_);

    /*
     * is_finished
     *
     * DESCRIPTION: make private finished-flag accessible
     * 
     * RETURNS: finished (which is true after first call to finish)
     */
    bool is_finished() {
        return finished;
    }

 private:
    /* internal functions */
    void initialise();
    void process_block(const unsigned char*);
    void get_result(void*);

    unsigned int A;                             /* accumulator 1 */
    unsigned int B;                             /* accumulator 2 */
    unsigned int C;                             /* accumulator 3 */
    unsigned int D;                             /* accumulator 4 */

    unsigned int message_length[2];             /* length of data */
    unsigned int stored_size;                   /* length of stored bytes */
    unsigned char stored[md5::BLOCK_SIZE * 2];  /* stored bytes */

    bool finished;                              /* object state */

    char signature[MD5_SIZE];                   /* stored signature */
    char str[MD5_STRING_SIZE];                  /* stored plain text hash */
};

/*
 * sig_to_string
 *
 * DESCRIPTION:
 * Convert a MD5 signature in a 16 byte buffer into a hexadecimal string representation.
 *
 * ARGUMENTS:
 * signature - a 16 byte buffer that contains the MD5 signature.
 * str - a string of characters which should be at least 33 bytes long (2 characters per MD5 byte and 1 for the \0).
 * str_len - the length of the string.
 */
extern void sig_to_string(const void* signature, char* str, const int str_len);

/*
 * sig_from_string
 *
 * DESCRIPTION:
 * Convert a MD5 signature from a hexadecimal string representation into a 16 byte buffer.
 *
 * ARGUMENTS:
 * signature - A 16 byte buffer that will contain the MD5 signature.
 * str - A string of charactes which _must_ be at least 32 bytes long (2 characters per MD5 byte).
 */
extern void sig_from_string(void* signature, const char* str);

}  // namespace md5



// Markus Iser: The following class just wraps all necessary stuff from above to produce md5 hashes from sequences of characters

#include <string>

class MD5 {
    md5::md5_t hasher;

 public:
    MD5() : hasher() { }
    ~MD5() { }

    void consume(const char* str, unsigned length) {
        // std::cout << std::string(str, length) << std::endl;
        hasher.process(str, length);
    }

    std::string produce() {
        unsigned char sig[MD5_SIZE];
        char str[MD5_STRING_SIZE];
        hasher.finish(sig);
        md5::sig_to_string(sig, str, sizeof(str));
        return std::string(str);
    }
};


#endif  // LIB_MD5_MD5_H_

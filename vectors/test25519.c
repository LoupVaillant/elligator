// This file is dual-licensed.  Choose whichever licence you want from
// the two licences listed below.
//
// The first licence is a regular 2-clause BSD licence.  The second licence
// is the CC-0 from Creative Commons. It is intended to release Monocypher
// to the public domain.  The BSD licence serves as a fallback option.
//
// SPDX-License-Identifier: BSD-2-Clause OR CC0-1.0
//
// ------------------------------------------------------------------------
//
// Copyright (c) 2022, Loup Vaillant
// All rights reserved.
//
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the
//    distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// ------------------------------------------------------------------------
//
// Written in 2022 by Loup Vaillant
//
// To the extent possible under law, the author(s) have dedicated all copyright
// and related neighboring rights to this software to the public domain
// worldwide.  This software is distributed without any warranty.
//
// You should have received a copy of the CC0 Public Domain Dedication along
// with this software.  If not, see
// <https://creativecommons.org/publicdomain/zero/1.0/>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "monocypher.h"

typedef int8_t   i8;
typedef uint8_t  u8;
typedef uint32_t u32;
typedef int32_t  i32;
typedef int64_t  i64;
typedef uint64_t u64;

#define FOR(i, start, end) for (size_t i = (start); i < (end); i++)

static int to_digit(int c)
{
    return c >= '0' && c <= '9' ? c - '0'
         : c >= 'a' && c <= 'f' ? c - 'a' + 10
         : c >= 'A' && c <= 'F' ? c - 'A' + 10
         : -1;
}

static int is_digit(int c)
{
    return to_digit(c) != -1;
}

int read_vector(u8 *vec, FILE *f)
{
    char c = getc(f);
    while (!is_digit(c) && c != ':' && c != EOF) {
        c = getc(f);
    }
    if (c == EOF ) {
        return 0;
    }
    while (c != ':' && c != EOF) {
        uint8_t msb = to_digit(c);  c = getc(f);
        uint8_t lsb = to_digit(c);  c = getc(f);
        *vec = lsb | (msb << 4);
        vec++;
    }
    if (c != ':') {
        fprintf(stderr, "Read_vector() error. Expected ':'\n");
        exit(1);
    }
    return 1;
}

void conclude(int status, char *test, unsigned nb_tests)
{
    fprintf(stderr, "%s (%3u tests): %s\n",
            status ? "FAILED" : "OK", nb_tests, test);
}

int main(int argc, const char *argv[])
{
    if (argc < 4) {
        fprintf(stderr, "Usage:\n  "
                "test25519.out direct.vec reverse.vec scalarmult.vec\n");
        exit(-1);
    }
    {
        FILE *f = fopen(argv[1], "r");
        u8 r[32];
        u8 u[32];
        u8 v[32]; // ignored
        u8 u2[32];
        int status = 0;
        unsigned nb_tests = 0;
        while (read_vector(r, f)) {
            read_vector(u, f);
            read_vector(v, f); // ignored
            crypto_hidden_to_curve(u2, r);
            status |= memcmp(u, u2, 32);
            nb_tests++;
        }
        conclude(status, "direct map", nb_tests);
        fclose(f);
    }
    {
        FILE *f = fopen(argv[2], "r");
        u8 u[32];
        u8 n[32];
        u8 s[32];
        u8 r[32];
        u8 r2[32];
        int status = 0;
        unsigned nb_tests = 0;
        while (read_vector(u, f)) {
            read_vector(n, f);
            read_vector(s, f);
            read_vector(r, f);
            int failed = crypto_curve_to_hidden(r2, u, n[0]);
            status |= (u8)failed != s[0];
            if (!failed) {
                status |= memcmp(r, r2, 32);
            }
            nb_tests++;
        }
        conclude(status, "reverse map", nb_tests);
        fclose(f);
    }
    {
        FILE *f = fopen(argv[3], "r");
        u8 sk [32];
        u8 pk [32];
        u8 pk2[32];
        int status = 0;
        unsigned nb_tests = 0;
        while (read_vector(sk, f)) {
            read_vector(pk, f);
            crypto_x25519_dirty_fast(pk2, sk);
            status |= memcmp(pk, pk2, 32);
            nb_tests++;
        }
        conclude(status, "scalarmult", nb_tests);
        fclose(f);
    }

    return 0;
}

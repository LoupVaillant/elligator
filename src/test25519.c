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
        u8 n[32]; // ignored
        u8 u2[32];
        int status = 0;
        unsigned nb_tests = 0;
        while (read_vector(r, f)) {
            read_vector(u, f);
            read_vector(v, f); // ignored
            read_vector(n, f); // ignored
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
        u8 p[32];
        u8 s[32];
        u8 r[32];
        u8 r2[32];
        int status = 0;
        unsigned nb_tests = 0;
        while (read_vector(u, f)) {
            read_vector(n, f);
            read_vector(p, f);
            read_vector(s, f);
            read_vector(r, f);
            int success = crypto_curve_to_hidden(r2, u, n[0] | (p[0] << 6));
            status |= (u8)success != s[0];
            status |= memcmp(r, r2, 32);
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

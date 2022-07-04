# This file is dual-licensed.  Choose whichever licence you want from
# the two licences listed below.
#
# The first licence is a regular 2-clause BSD licence.  The second licence
# is the CC-0 from Creative Commons. It is intended to release Monocypher
# to the public domain.  The BSD licence serves as a fallback option.
#
# SPDX-License-Identifier: BSD-2-Clause OR CC0-1.0
#
# ------------------------------------------------------------------------
#
# Copyright (c) 2022, Loup Vaillant
# All rights reserved.
#
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the
#    distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ------------------------------------------------------------------------
#
# Written in 2022 by Loup Vaillant
#
# To the extent possible under law, the author(s) have dedicated all copyright
# and related neighboring rights to this software to the public domain
# worldwide.  This software is distributed without any warranty.
#
# You should have received a copy of the CC0 Public Domain Dedication along
# with this software.  If not, see
# <https://creativecommons.org/publicdomain/zero/1.0/>

CC            = gcc -std=c99
CFLAGS        = -Wall -Wextra -pedantic -O3
PYFILES       = core.py elligator.py gen_vectors.py
VECTORS_25519 = ../vectors/curve25519_direct.vec  \
                ../vectors/curve25519_inverse.vec \
                ../vectors/curve25519_scalarmult.vec
VECTORS_448   = ../vectors/curve448_direct.vec  \
                ../vectors/curve448_inverse.vec \
                ../vectors/curve448_scalarmult.vec

.PHONY: all test clean

all: $(VECTORS_25519) $(VECTORS_448) ../vectors/hash_to_curve25519.vec

test: test25519 test_hash_to_25519

test25519: test25519.out $(VECTORS_25519)
	./$^

test_hash_to_25519: test_hash_to_25519.out ../vectors/hash_to_curve25519.vec
	./$^

../vectors/curve25519_%.vec: $(PYFILES) curve25519.py
	./gen_vectors.py curve25519 $* >$@

../vectors/curve448_%.vec: $(PYFILES) curve448.py
	./gen_vectors.py curve448 $* >$@

../vectors/hash_to_curve25519.vec: hash_to_curve25519.py
	./$< >$@

test25519.out: test25519.c
	$(CC) $(CFLAGS) $< -o $@               \
            $$(pkg-config monocypher --cflags) \
            $$(pkg-config monocypher --libs)

test_hash_to_25519.out: test_hash_to_25519.c
	$(CC) $(CFLAGS) $< -o $@              \
            $$(pkg-config libsodium --cflags) \
            $$(pkg-config libsodium --libs)

clean:
	rm -f *.out ../vectors/*.vec

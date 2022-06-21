#! /usr/bin/env python3

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

from curve25519 import *
from elligator import *
from random    import randrange
from random    import seed
import hashlib

def map_to_curve(random):
    """Maps a uniform random 256 bit number to a curve point

    This is compatible with libsodium.

    Contrary to what the Hash to Curve RFC recommends, we don't use a
    big hash to reduce it modulo p.  Instead we just chop off one
    bit. The resulting field is close enough to a power of two that the
    deviation from perfect randomness is undetectable.
    """
    y_sign  = random // 2**255            # Get sign of Edwards x coordinate
    r       = random %  2**255            # Elligator representative
    u, _    = dir_map(GF(r))              # Ignore Montgomery v coordinate
    x, y, z = to_edwards(u)               # Convert to Edwards
    if x.to_num() % 2 != y_sign: x = -x   # Set sign of Edwards x coordinate
    x, y, z = Ed.scalarmult((x, y, z), 8) # Multiply by cofactor

    # Serialise Edwards point (divide, get sign of x)
    z       = z.invert()
    y       = y * z
    x       = x * z
    x_sign  = x.to_num() % 2              # Negative means odd here.
    point   = y.to_num() + x_sign * 2**255
    return point

# Generate the actual test vectors, print them in stdout.
seed(12345)  # cheap determinism for the random test vectors
vectors = []
for i in range(64):
    r = randrange(2**256)
    p = map_to_curve(r)
    vectors.append(vectors_to_string([r, p]))
print("\n\n".join(vectors))

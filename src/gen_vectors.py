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
# Copyright (c) 2020-2022, Loup Vaillant
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
# Written in 2020-2022 by Loup Vaillant
#
# To the extent possible under law, the author(s) have dedicated all copyright
# and related neighboring rights to this software to the public domain
# worldwide.  This software is distributed without any warranty.
#
# You should have received a copy of the CC0 Public Domain Dedication along
# with this software.  If not, see
# <https://creativecommons.org/publicdomain/zero/1.0/>

import sys

# collect arguments
if len(sys.argv) < 3:
    raise ValueError('Usage: gen_vectors.py curve vectors')
curve   = sys.argv[1]
vectors = sys.argv[2]

# Import curve module
if curve == "curve25519": from curve25519 import *
else: raise ValueError('')

# remaining imports
from elligator import *
from random    import randrange
from random    import seed


##############
# Direct map #
##############
def direct_map_vectors(r, pad):
    """One set of test vectors for the direct map"""
    p1   = map(r)
    p2   = map(-r)
    p3   = fast_map(r)
    u, v = p1
    r1   = rev_map     (u, v.is_negative())
    r2   = fast_rev_map(u, v.is_negative())
    if p1 != p2: raise ValueError('p mismatch (map vs fast_map)')
    if p1 != p3: raise ValueError('p mismatch (r vs -r)')
    if r1 != r2: raise ValueError('r mismatch (rev_map vs fast_rev_map')
    if r  != r1: raise ValueError("r mismatch (round trip failure)")

    return vectors_to_string([
        r.to_num() + pad * 2**GF.msb,
        u,
        v,
        v.is_negative()
    ])

def direct_map_all_vectors():
    """All test vectors for the direct map"""
    seed(12345)  # cheap determinism for the random test vectors
    vectors = []
    max_pad = GF.max_pad * 2 # because all representatives are positive

    # Representative 0 maps to the point 0, 0
    for pad in range(max_pad):
        vectors.append(direct_map_vectors(GF(0), pad))

    # Random representatives map to their respective point
    for i in range(256):
        pad = i % max_pad
        r   = GF(randrange(0, GF.p - 1)).abs()
        vectors.append(direct_map_vectors(r, pad))

    return "\n\n".join(vectors)


###############
# Reverse map #
###############
def reverse_map(u, v_is_negative):
    r1 = rev_map     (u, v_is_negative)
    r2 = fast_rev_map(u, v_is_negative)
    if r1 != r2: raise ValueError('r mismatch (rev_map vs fast_rev_map)')
    return r1

def on_curve():
    u = GF(randrange(0, GF.p - 1))
    while not is_square(u**3 + A * u**2 + B * u):
        u = GF(randrange(0, GF.p - 1))
    return u

def reverse_map_vectors_fail(pad):
    """Test vectors for 2 opposite failing points"""
    p = "%02x:" % (pad)
    u = on_curve()
    r = reverse_map(u, False)
    while not r is None:
        u = on_curve()
        r = reverse_map(u, False)
    if not reverse_map(u, True) is None:
        raise ValueError('Reverse map should fail')
    v_pos = vectors_to_string([u, False, p, "ff:", ":"])
    v_neg = vectors_to_string([u, True , p, "ff:", ":"])
    return v_pos + "\n\n" + v_neg

def reverse_map_vectors_ok(pad):
    """Test vectors for 2 opposite succeeding points"""
    p = "%02x:" % (pad)
    u  = on_curve()
    rp = reverse_map(u, False)
    while rp is None:
        u  = on_curve()
        rp = reverse_map(u, False)
    rn = reverse_map(u, True)
    if rn is None: raise ValueError('Reverse map should succeed')
    v_pos = vectors_to_string([u, False, p, "00:", rp.to_num()+pad*2**GF.msb])
    v_neg = vectors_to_string([u, True , p, "00:", rn.to_num()+pad*2**GF.msb])
    return v_pos + "\n\n" + v_neg

def reverse_map_all_vectors():
    """All test vectors for the reverse map"""
    seed(12345)  # cheap determinism for the random test vectors
    vectors = []
    max_pad = GF.max_pad * 2 # because all representatives are positive

    # point (0, 0) maps to representative 0
    for pad in range(max_pad):
        p = "%02x:" % (pad)
        vectors.append(vectors_to_string([0, False, p,"00:", pad*2**GF.msb]))

    # some points that do not map
    for i in range(16):
        pad = i % max_pad
        vectors.append(reverse_map_vectors_fail(pad))

    # lots of points that do map
    for i in range(256):
        pad = i % max_pad
        vectors.append(reverse_map_vectors_ok(pad))

    return "\n\n".join(vectors)


##############
# Scalarmult #
##############
def scalarmult_vectors(scalar, c):
    scalar -= scalar % 8
    scalar += c
    return vectors_to_string([
        scalar,
        scalarmult(scalar, c)
    ])

def scalarmult_all_vectors():
    """All test vectors for scalar multiplication"""
    seed(12345)  # cheap determinism for the random test vectors
    vectors = []
    for i in range(64):
        scalar = randrange(2**(GF.nb_bytes * 8))
        c      = i % cofactor
        vectors.append(scalarmult_vectors(scalar, c))
    return "\n\n".join(vectors)


################
# Main program #
################
vectors_map = {"direct"    : direct_map_all_vectors,
               "reverse"   : reverse_map_all_vectors,
               "scalarmult": scalarmult_all_vectors,
               }
print(vectors_map[vectors]())

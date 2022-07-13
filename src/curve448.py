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

import core
from core import *

##################
# Isogeny switch #
##################

# Traditionally, Ed448 doesn't use X448's birationally equivalent curve,
# but an isogeny thereof.  This introduces a slight complication when
# trying to lean on existing Ed448 code.
#
# Here I decided to describe both alternatives:
# - If isogeny is True, we use the standard Ed448 code.
# - If isogeny is False, we use the birationally equivalent Edwards curve.
#
# Both methods generate the exact same results.
isogeny = True


####################
# field parameters #
####################
GF.p = 2**448 - 2**224 - 1

def is_negative(self):
    """True iff self is odd"""
    return (self.val % GF.p) % 2 == 1

GF.is_negative = is_negative


#########################
# Square root functions #
#########################
def sqrt(n):
    """Principal square root of n

    If n is not a square, the behaviour is undefined.

    To know how it works, refer to https//elligator.org/sqrt
    """
    square = n**((GF.p+1) // 4)
    square = cmove(square, -square, square.is_negative())
    return square

def inv_sqrt(x):
    """Inverse square root of x

    Returns (0         , True ) if x is zero.
    Returns (sqrt( 1/x), True ) if x is non-zero square.
    Returns (sqrt(-1/x), False) if x is not a square.

    The return value is *not* guaranteed to be non-negative.
    """
    isr       = x**((GF.p - 3) // 4)
    legendre  = x * isr**2         # x**((p-1)/2) == -1, 0, or 1
    is_square = legendre != GF(-1) # use constant time comparison
    return isr, is_square

core.sqrt     = sqrt
core.inv_sqrt = inv_sqrt


################################
# Birational map (and isogeny) #
################################
def mt_to_edwards(u):
    y = (u + GF(1)) / (u - GF(1))
    x = sqrt((y**2 - GF(1)) / (Ed.d * y**2 - GF(1)))
    return (x, y, GF(1))

def isogeny_to_ed(point):
    x, y, z = point
    x2 = x**2
    y2 = y**2
    du = z**2*GF(2) - x2 - y2
    v  = y2 - x2 # dv, actually. We use v to save space
    d  = du * v
    u  = v * x * y * GF(2)*sqrt(Ed.d)
    v  = (y2 + x2) * du
    return (u, v, d)

def edwards_to_mt(point):
    x, y, z = point  # in projective coordinates
    return (y + z) / (y - z)


####################
# Curve parameters #
####################

# Montgomery constants (We already assume B = 1)
Mt.A = GF(156326)

# Edwards constants
Ed.a = GF(1) # 1 -> not twisted
if isogeny: Ed.d = GF(-39081)
else      : Ed.d = GF(39082) / GF(39081)

# curve order and cofactor
Mt.order    = 2**446-0x8335dc163bb124b65129c96fde933d8d723a70aadc873d6d54a7bb0d
Mt.cofactor = 4

# Standard base point, that generates the prime order sub-group
Mt.base = GF(5)
if isogeny:
    # From RFC 8032
    Ed.base = (GF(224580040295924300187604334099896036246789641632564134246125461686950415467406032909029192869357953282578032075146446173674602635247710),
               GF(298819210078481492676017930443930673437544040154080242095928241372331506189835876003536878655418784733982303233503462500531545062832660),
               GF(1))
else:
    # Base point of the (non-standard) birational curve
    Ed.base = mt_to_edwards(Mt.base)

# Low order point (of order 4), used to add the cofactor component
# There are 2 such points: (1, 0) and (-1, 0)
# We chose (1, 0) somewhat arbitrarily
Ed.lop = (GF(1), GF(0), GF(1))

# "Dirty" Base point, that generates the whole curve.
# mt_base_c = mt_base + (lop * co_clear)
co_clear  = Mt.order % Mt.cofactor # 3
lop_c     = Ed.scalarmult(Ed.lop, co_clear)
ed_base   = isogeny_to_ed(Ed.base) if isogeny else Ed.base
Ed.base_c = Ed.add(ed_base, lop_c)
Mt.base_c = edwards_to_mt(Ed.base_c)

def add_lop(point, i):
    """Adding a low order point, fast

    # Equivalent to the following, except constant time
    x, y, z = point
    if i == 0: return  x,  y, z
    if i == 1: return  y, -x, z
    if i == 2: return -x, -y, z
    if i == 3: return -y,  x, z
    """
    x, y, z = point
    l    = (i//1) % 2 == 1
    h    = (i//2) % 2 == 1
    x, y = cswap(x, y, l)
    x    = cmove(x, -x, h)
    y    = cmove(y, -y, l != h) # use XOR instead of !=
    return x, y, z

# Replacing the generic Ed.co_scalarmult()
# with a custom method that uses the fast add_lop()
# instead of the regular select_lop() then add().
# This saves a full point addition.
def co_scalarmult(scalar, c):
    main_point  = Ed.scalarmult(Ed.base, clamp(scalar))
    if isogeny:
        main_point = isogeny_to_ed(main_point)
    low_order_p = Ed.scalarmult(Ed.lop, c)
    montgomery1 = edwards_to_mt(Ed.add(main_point, low_order_p)) # slow
    montgomery2 = edwards_to_mt(add_lop(main_point, c))          # fast
    if montgomery1 != montgomery2:
        raise ValueError('Incoherent low order point selection')
    return montgomery2

Ed.co_scalarmult = co_scalarmult


########################
# Elligator parameters #
########################
core.Z       = GF(-1)
core.ufactor = -core.Z            # ufactor ==  1
core.vfactor = sqrt(core.ufactor) # vfactor == -1

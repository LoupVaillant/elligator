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

import core
from core import *

####################
# field parameters #
####################
GF.p = 2**255 - 19

def is_negative(self):
    """True iff self is in [p.+1 / 2.. p-1]

    An alternative definition is to just test whether self is odd.
    """
    dbl = (self.val * 2) % GF.p
    return dbl % 2 == 1

GF.is_negative = is_negative


#########################
# Square root functions #
#########################
sqrt_m1 = (GF(2)**((GF.p-1) // 4)).abs() # sqrt(-1)

def sqrt(n):
    """ Non-negative square root of n

    If n is not a square, the behaviour is undefined.

    To know how it works, refer to https//elligator.org/sqrt
    """
    root = n**((GF.p+3) // 8)
    root = cmove(root, root * sqrt_m1, root * root != n)
    return root.abs() # returns the non-negative square root

def inv_sqrt(x):
    """Inverse square root of x

    Returns (0               , True ) if x is zero.
    Returns (sqrt(1/x)       , True ) if x is non-zero square.
    Returns (sqrt(sqrt(-1)/x), False) if x is not a square.

    The return value is *not* guaranteed to be non-negative.
    """
    isr       = x**((GF.p - 5) // 8)
    quartic   = x * isr**2
    # Use constant time comparisons in production code
    m_sqrt_m1 = quartic == GF(-1) or quartic == -sqrt_m1
    is_square = quartic == GF(-1) or quartic == GF(1) or x == GF(0)
    isr       = cmove(isr, isr * sqrt_m1, m_sqrt_m1)
    return isr, is_square

core.sqrt     = sqrt
core.inv_sqrt = inv_sqrt


##################
# Birational map #
##################
def to_edwards(u):
    y = (u - GF(1)) / (u + GF(1))
    x = sqrt((y**2 - GF(1)) / (Ed.d * y**2 + GF(1)))
    return (x, y, GF(1))

def to_montgomery(point):
    x, y, z = point  # in projective coordinates
    return (z + y) / (z - y)

Ed.to_mt = to_montgomery


####################
# Curve parameters #
####################

# Montgomery constants (We already assume B = 1)
Mt.A = GF(486662)

# Twisted Edwards constants
Ed.a = GF(-1)
Ed.d = GF(-121665) / GF(121666)

# curve order and cofactor
Mt.order    = 2**252 + 27742317777372353535851937790883648493
Mt.cofactor = 8

# Standard base point, that generates the prime order sub-group
Mt.base = GF(9)                # Montgomery base point
Ed.base = to_edwards(Mt.base)  # Edwards base point

# Low order point (of order 8), used to add the cofactor component
# There are 4 such points, that differ only by the sign of
# their coordinates: (x, y), (x, -y), (-x, y), (-x, -y)
# We chose the one whose both coordinates are positive (below GF.p // 2)
lop_x       = sqrt((sqrt(Ed.d + GF(1)) + GF(1)) / Ed.d)
lop_y       = -lop_x * sqrt_m1
Ed.lop = (lop_x, lop_y, GF(1))

# "Dirty" Base point, that generates the whole curve.
# Mt.base_c = Mt.base + (lop * co_clear)
co_clear  = Mt.order % Mt.cofactor # 5
lop_c     = Ed.scalarmult(Ed.lop, co_clear)
Ed.base_c = Ed.add(Ed.base, lop_c)
Mt.base_c = Ed.to_mt(Ed.base_c)

# Constant time selection of the low order point
# Using tricks to minimise the size of the look up table
# It's a faster alternative to scalar multiplication.
#
# The 8 low order points are as follows:
#
# [0]lop = ( 0       ,  1    )
# [1]lop = ( lop_x   ,  lop_y)
# [2]lop = ( sqrt(-1), -0    )
# [3]lop = ( lop_x   , -lop_y)
# [4]lop = (-0       , -1    )
# [5]lop = (-lop_x   , -lop_y)
# [6]lop = (-sqrt(-1),  0    )
# [7]lop = (-lop_x   ,  lop_y)
#
# The select() subroutine takes advantage of the common cyclical
# patterns in the values of the x and y coordinates.
def select_lop(i):
    def select(x, k, i):
        r = GF(0)
        r = cmove(r,  k, (i // 2) % 2 == 1) # bit 1
        r = cmove(r,  x, (i // 1) % 2 == 1) # bit 0
        r = cmove(r, -r, (i // 4) % 2 == 1) # bit 2
        return r
    x = select(lop_x, sqrt_m1, i  )
    y = select(lop_y, GF(1)  , i+2)
    return (x, y, GF(1))

Ed.select_lop = select_lop


########################
# Elligator parameters #
########################
core.Z       = GF(2)               # sqrt(-1) is sometimes faster...
core.ufactor = -core.Z * sqrt_m1   # ...because then both ufactor
core.vfactor = sqrt(core.ufactor)  # and vfactor are equal to 1

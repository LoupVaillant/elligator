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
    chi       = x * isr**2    # x**((p-1)/2) == -1, 0, or 1
    is_square = chi != GF(-1) # use constant time comparison
    return isr, is_square

core.sqrt     = sqrt
core.inv_sqrt = inv_sqrt


#################################################
# Birational map between Edwards and Montgomery #
#################################################
def from_edwards(point):
    x, y, z = point  # in projective coordinates
    return (y + z) / (y - z)

def to_edwards(u):
    y = (u + GF(1)) / (u - GF(1))
    x = sqrt((y**2 - GF(1)) / (core.D_e * y**2 - GF(1)))
    return (x, y, GF(1))

core.from_edwards = from_edwards
core.to_edwards   = to_edwards


####################
# curve parameters #
####################
# Montgomery constants
core.A = GF(156326)
core.B = GF(1)

# Edwards constants
core.A_e = GF(1) # 1 -> not twisted
core.D_e = GF(39082) / GF(39081)

# curve order and cofactor
# co_clear is choosen such that for all x:
# - (x * order * co_clear) % order    = 0
# - (x * order * co_clear) % cofactor = x % cofactor
# The goal is to preserve the cofactor and eliminate the main factor.
core.order    = 2**446-0x8335dc163bb124b65129c96fde933d8d723a70aadc873d6d54a7bb0d
core.cofactor = 4
core.co_clear = 3

# Low order point (of order 4), used to add the cofactor component
# There are 2 such points: (1, 0) and (-1, 0)
# We chose the one with positive coordinates.
core.lop = (GF(1), GF(0), GF(1))

# Standard base point, that generates the prime order sub-group
core.mt_base = GF(5)                     # Montgomery base point
core.ed_base = to_edwards(core.mt_base)  # Edwards base point

# "Dirty" Base point, that generates the whole curve.
# mt_base_c = mt_base + (lop * co_clear)
lop_c          = ed_scalarmult(core.lop, core.co_clear)
core.ed_base_c = point_add(core.ed_base, lop_c)
core.mt_base_c = core.from_edwards(core.ed_base_c)

# Constant time selection of the low order point
# Using tricks to minimise the size of the look up table
# It's a faster alternative to scalar multiplication.
#
# The 4 low order points are as follows:
#
# [0]lop = ( 0,  1)
# [1]lop = ( 1, -0)
# [2]lop = (-0, -1)
# [3]lop = (-1,  0)
#
# The select() subroutine takes advantage of the common cyclical
# patterns in the values of the x and y coordinates.
def select_lop(i):
    def select(i):
        r = GF(0)
        r = cmove(r, GF(1), (i // 1) % 2 == 1) # bit 0
        r = cmove(r, -r   , (i // 2) % 2 == 1) # bit 1
        return r
    x = select(i  )
    y = select(i+1)
    return (x, y, GF(1))

core.select_lop = select_lop

########################
# Elligator parameters #
########################
core.Z       = GF(-1)
core.ufactor = -core.Z            # ufactor ==  1
core.vfactor = sqrt(core.ufactor) # vfactor == -1

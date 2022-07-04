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
# Copyright (c) 2020-2021, Loup Vaillant
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
# Written in 2020-2021 by Loup Vaillant
#
# To the extent possible under law, the author(s) have dedicated all copyright
# and related neighboring rights to this software to the public domain
# worldwide.  This software is distributed without any warranty.
#
# You should have received a copy of the CC0 Public Domain Dedication along
# with this software.  If not, see
# <https://creativecommons.org/publicdomain/zero/1.0/>

import math # log

####################
# Field arithmetic #
####################
class GF():
    """Finite field over some prime number.

    Only prime fields are supported.
    No extension field here.
    This class is supposed to be derived,
    so the prime number p is defined

    the fowlowing is not implemented, and must be defined
    with inheritance or monkey patching:
    - p                 : characteristic of the field
    - is_negative(self) : set of negative field elements
    """

    def __init__(self, x):
        GF.msb         = round(math.log(GF.p, 2)) - 1
        GF.nb_bytes    = math.ceil((GF.msb + 1) / 8)
        GF.nb_pad_bits = GF.nb_bytes * 8 - GF.msb - 1
        GF.max_pad     = 2**GF.nb_pad_bits
        self.val       = x % self.p

    # Basic arithmetic operations
    def __neg__     (self   ): return GF(-self.val                            )
    def __add__     (self, o): return GF( self.val +  o.val                   )
    def __sub__     (self, o): return GF( self.val -  o.val                   )
    def __mul__     (self, o): return GF((self.val *  o.val         ) % self.p)
    def __truediv__ (self, o): return GF((self.val *  o.invert().val) % self.p)
    def __floordiv__(self, o): return GF( self.val // o) # same as __truediv__
    def __pow__     (self, s): return GF(pow(self.val, s       , self.p))
    def invert      (self   ): return GF(pow(self.val, self.p-2, self.p))

    def __eq__(self, other): return self.val % self.p == other.val % self.p
    def __ne__(self, other): return self.val % self.p != other.val % self.p

    def is_positive(self)  : return not self.is_negative()
    def abs(self):
        if self.is_positive(): return  self
        else                 : return -self

    def to_num(self):
        return self.val % self.p

    def __str__ (self): return str(self.to_num())
    def __repr__(self): return str(self.to_num())

def to_hex(n):
    """Converts a number in hexadecimal (little endian)"""
    p = GF.p
    s = ""
    while p != 0:
        s += format(n % 256, '02x')
        n //= 256
        p //= 256
    return s + ':'

def vectors_to_string(values):
    strings = []
    for v in values:
        if   type(v) is GF  : strings.append(to_hex(v.to_num()))
        elif type(v) is bool: strings.append("01:" if v else "00:")
        elif type(v) is str : strings.append(v)
        else                : strings.append(to_hex(v))
    return "\n".join(strings)


##################################
# Scalar clamping (X25519, X448) #
##################################
def clamp(scalar):
    clamped = scalar - scalar % Mt.cofactor
    clamped = clamped % 2**GF.msb
    clamped = clamped + 2**GF.msb
    return clamped


################################
# Basic square root operations #
################################
def legendre(n):
    """Legendre symbol:

    returns  0 if n is zero
    returns  1 if n is a non-zero square
    returns -1 if n is not a square
    """
    return n**((GF.p-1)//2)

def is_square(n):
    c = legendre(n)
    return c == GF(0) or c == GF(1)


###########################
# Constant time selection #
###########################
def cswap(a, b, swap):
    """Conditionnal swapping

    Usage:
        a, b = cswap(a, b, swap)

    Swaps a and b if swap is true,
    does nothing otherwise.

    Production code is supposed to run in constant time.
    Be aware that this code does not.
    """
    if swap: return b, a
    else   : return a, b

def cmove(a, b, move):
    """Conditionnal assignment

    Usage:
        a = cmove(a, b, move)

    Assigns the value of b to a if move is true,
    does nothing otherwise.

    Production code is supposed to run in constant time.
    Be aware that this code does not.
    """
    if move: return b
    else   : return a


#################
# Edwards curve #
#################
class Ed():
    """Edwards curve (possibly twisted)

    Not really a class.  Think of it as a namespace.
    The following must be defined with monkey patching:
    - a         : curve constant
    - d         : curve constant
    - lop       : low order point
    - to_mt     : convertion function from Edwards to montgomery
    - select_lop: fast low order point selection
    """
    def add(p1, p2):
        """Point addition, using projective coordinates

        Formula with affine coordinates:
            denum = d*x1*x2*y1*y2
            x     = (x1*y2 +   x2*y1) / (1 + denum)
            y     = (y1*y2 - a*x1*x2) / (1 - denum)

        We use projective coordinates to avoid expensive divisions:
            P = (X, Y, Z)
            x = X / Z
            y = Y / Z
        """
        x1, y1, z1 = p1
        x2, y2, z2 = p2
        denum = Ed.d*x1*x2*y1*y2
        t1    = z1 * z2
        t2    = t1**2
        xt    = t1 * (x1*y2 +        x2*y1)
        yt    = t1 * (y1*y2 - Ed.a*x1*x2)
        zx    = t2 + denum
        zy    = t2 - denum
        return (xt*zy, yt*zx, zx*zy)

    def check_point(p):
        """Is the point on the curve?

        Check that the point p (in projective coordinates to avoid
        expensive divisions) matches the following twisted Edwards
        equation: a*x^2 + y^2 = 1 + d*x^2*y^2
        """
        x, y, z = p
        x2, y2, z2, z4 = (x**2, y**2, z**2, z**4)
        if (Ed.a*x2 + y2)*z2 != z4 + Ed.d*x2*y2:
            raise ValueError("Point not on the curve!!")

    def scalarmult(point, scalar):
        """Scalar multiplication in Edwards space"""
        Ed.check_point(point)
        acc    = (GF(0), GF(1), GF(1))
        binary = [int(c) for c in list(format(scalar, 'b'))]
        for i in binary:
            acc = Ed.add(acc, acc)
            Ed.check_point(acc)
            if i == 1:
                acc = Ed.add(acc, point)
                Ed.check_point(acc)
        return acc

    def co_scalarmult(scalar, c):
        """Scalarmult with cofactor

        Returns a point converted to Montgomery.
        There are two equivalent ways to select the low order point:
        - Scalar multiplication
        - Constant time selection from a table
        Selecting from a table is generally very simple and fast
        """
        main_point  = Ed.scalarmult(Ed.base, clamp(scalar))
        low_order1  = Ed.scalarmult(Ed.lop, c)
        low_order2  = Ed.select_lop(c)
        montgomery1 = Ed.to_mt(Ed.add(main_point, low_order1))
        montgomery2 = Ed.to_mt(Ed.add(main_point, low_order2))
        if montgomery1 != montgomery2:
            raise ValueError('Incoherent low order point selection')
        return montgomery1


####################
# Montgomery curve #
####################
class Mt():
    """Montgomery curve

    The following must be defined with monkey patching:
    - A     : curve constant
    - base_c: special base point that covers the whole curve

    The curve constant B is assumed equal to 1 (it has to be for the
    Montgomery curve to be compatible with Elligator2).
    """
    def scalarmult(u, scalar):
        """Scalar multiplication in Montgomery space

        This is an "X-only" laddder, that only uses the u coordinate.
        This conflates points (u, v) and (u, -v).
        """
        u2, z2 = GF(1), GF(0) # "zero" point
        u3, z3 = u    , GF(1) # "one"  point
        binary = [int(c) for c in list(format(scalar, 'b'))]
        for b in binary:
            # Montgomery ladder step:
            # if b == 0, then (P2, P3) == (P2*2 , P2+P3)
            # if b == 1, then (P2, P3) == (P2+P3, P3*2 )
            swap   = b == 1  # Use constant time comparison
            u2, u3 = cswap(u2, u3, swap)
            z2, z3 = cswap(z2, z3, swap)
            u3, z3 = ((u2*u3 - z2*z3)**2, (u2*z3 - z2*u3)**2 * u)
            u2, z2 = ((u2**2 - z2**2)**2,
                      GF(4)*u2*z2*(u2**2 + Mt.A*u2*z2 + z2**2))
            u2, u3 = cswap(u2, u3, swap)
            z2, z3 = cswap(z2, z3, swap)
        return u2 / z2

    def co_scalarmult(scalar, c):
        """Scalarmult with cofactor"""
        co_cleared = (c % Mt.cofactor) * Mt.order  # cleared main factor
        combined   = clamp(scalar) + co_cleared
        return Mt.scalarmult(Mt.base_c, combined)


############################
# Scalarmult with cofactor #
############################

# Keeping a random cofactor is important to keep points
# indistinguishable from random.  (Else we'd notice all representatives
# represent points with cleared cofactor.  Not exactly random.)

# Perform the different scalar multiplications and compare them.
# All methods are supposed to yield the same results.
def co_scalarmult(scalar, c):
    p1 = Ed.co_scalarmult(scalar, c)
    p2 = Mt.co_scalarmult(scalar, c)
    if p1 != p2:
        raise ValueError('Incoherent scalarmult')
    return p1

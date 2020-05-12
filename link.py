#!/usr/bin/env python3
# -*- coding: utf-8 -*-


class Link:
    def __init__(self, seg1, seg2, hicVal=0):
        self.ends = [seg1, seg2]
        self.hicValue = hicVal

# getters

def _get_ends(self):
    return self._ends


def _get_hicValue(self):
    return self._hicValue


# setters


def _set_ends(self, seg1, seg2):
    self._ends = [seg1, seg2]


def _set_hicValue(self, hicVal):
    self._hicValue = hicVal


# properties

ends = property(_get_ends, _set_ends)
hicValue = property(_get_hicValue, _set_hicValue)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-


class Segment:
    def __init__(self, segId, segSeq):
        self.id = segId
        self.seq = segSeq
        self.head = []
        self.tail = []


# getters


def _get_id(self):
    return self._id


def _get_seq(self):
    return self._seq


def _get_head(self):
    return self._head


def _get_tail(self):
    return self._tail


# setters


def _set_id(self, newId):
    self._id = str(newId)


def _set_seq(self, newSeq):
    self._seq = str(newSeq)


# properties

id = property(_get_id, _set_id)
seq = property(_get_seq, _set_seq)
head = property(_get_head)
tail = property(_get_tail)

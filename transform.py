#!/usr/bin/env python

import sys
from pycparser import parse_file, c_generator
from pycparser.c_ast import NodeVisitor
import pycparser.c_ast as c_ast
from itertools import zip_longest
import re

class BetterNodeVisitor(object):
    _method_cache = None

    def visit(self, node, parent, name, index):
        """ Visit a node.
        """

        if self._method_cache is None:
            self._method_cache = {}

        visitor = self._method_cache.get(node.__class__.__name__, None)
        if visitor is None:
            method = 'visit_' + node.__class__.__name__
            visitor = getattr(self, method, self.generic_visit)
            self._method_cache[node.__class__.__name__] = visitor

        return visitor(node, parent, name, index)

    def generic_visit(self, node, parent, name, index):
        """ Called if no explicit visitor function exists for a
            node. Implements preorder visiting of the node.
        """
        for i, (name, c) in enumerate(node.children()):
            self.visit(c, node, name, i)

def grouper(iterable, n, *, incomplete='fill', fillvalue=None):
    "Collect data into non-overlapping fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, fillvalue='x') --> ABC DEF Gxx
    # grouper('ABCDEFG', 3, incomplete='strict') --> ABC DEF ValueError
    # grouper('ABCDEFG', 3, incomplete='ignore') --> ABC DEF
    args = [iter(iterable)] * n
    if incomplete == 'fill':
        return zip_longest(*args, fillvalue=fillvalue)
    if incomplete == 'strict':
        return zip(*args, strict=True)
    if incomplete == 'ignore':
        return zip(*args)
    else:
        raise ValueError('Expected fill, strict, or ignore')

class ScopFinder(BetterNodeVisitor):
    def __init__(self):
        self.pragmas = []
        self.parents = []
        self.scops = []

    def visit_Pragma(self, node, parent, name, index):
        if "endscop" in node.string:
            self.pragmas.append((node, index))
        elif "scop" in node.string:
            self.pragmas.append((node, index))
            self.parents.append(parent)

    def get_scops(self):
        for i, (p1, p2) in enumerate(grouper(self.pragmas, 2)):
            parent = self.parents[i]
            children = list(iter(parent))
            scop = children[p1[1]:p2[1]+1]
            self.scops.append(scop)

        return self.scops

def replace_code_between_lines(indexed_lines, replacement, l0, l1):
    for i, (line_no, line) in enumerate(indexed_lines):
        if line_no == l0:
            indexed_lines[i] = (l0, replacement)
        elif line_no > l0 and line_no <= l1:
            indexed_lines[i] = (l0, "")

    return indexed_lines

def transform(source, destination):
    ast = parse_file(source, use_cpp=True, cpp_args=[r'-I../pycparser/utils/fake_libc_include',
                                                     r'-I.'])

    sf = ScopFinder()
    sf.visit(ast, None, None, None)
    scops = sf.get_scops()

    with open(source, "r") as s, open(destination, "w") as d:
        text = s.read()
        indexed_lines = list(enumerate(text.split("\n")))

        for s in scops:
            l0 = s[0].coord.line - 1
            l1 = s[-1].coord.line - 1
            generator = c_generator.CGenerator()
            replacement = generator.visit(c_ast.Compound(s))
            replace_code_between_lines(indexed_lines, replacement, l0, l1)

        lines = [il[1] for il in indexed_lines]
        d.write("\n".join(lines))

if __name__ == "__main__":
    if len(sys.argv) > 2:
        transform(sys.argv[1], sys.argv[2])
    else:
        print("Usage: transform.py source dest")

#!/usr/bin/env python

import sys
from pycparser import parse_file, c_generator
from pycparser.c_ast import NodeVisitor
import pycparser.c_ast as c_ast
from itertools import zip_longest
from utils import *
import re

source = "./gemver-pre.c"
destination = "./gemver.c"

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

        a_for_node = scops[0][1]

        lines = [il[1] for il in indexed_lines]
        d.write("\n".join(lines))

if __name__ == "__main__":
    if len(sys.argv) > 2:
        transform(sys.argv[1], sys.argv[2])
    else:
        print("Usage: transform.py source dest")

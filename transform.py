#!/usr/bin/env python

import sys
from pycparser import parse_file, c_generator
from pycparser.c_ast import NodeVisitor
import pycparser.c_ast as c_ast
from itertools import zip_longest
from utils import *
import re

class TransformationApplier(BetterNodeVisitor):
    def visit_Pragma(self, node, parent, name, index):
        if "unroll" in node.string:
            result = re.findall(r"([_a-zA-Z][_a-zA-Z0-9]*):(\d+)", node.string)
            unroll_guide = dict(result)
            unroll_guide = {k: int(v) for k, v in unroll_guide.items()}

            original_loop = parent.block_items[index + 1]
            print("ORIGINAL LOOP:")
            dump(original_loop)
            loop = preform_all_unrolling(original_loop, unroll_guide)
            if "licm" in node.string:
                loop = preform_all_hoisting(loop)
            if "prefetch" in node.string:
                loop = preform_all_prefetching(loop)
            print("NEW LOOP:")
            dump(loop)
            parent.block_items[index + 1] = loop

# source = "./benchmarks/fdtd-apml.ppcg.c.in"
# destination = "./nothing.txt"

def transform(source, destination):
    ast = parse_file(source, use_cpp=True, cpp_args=[r'-I./pycparser/utils/fake_libc_include',
                                                     r'-I.'])

    sf = ScopFinder()
    sf.visit(ast, None, None, None)
    scops = sf.get_scops()

    with open(source, "r") as s, open(destination, "w") as d:
        text = s.read()
        indexed_lines = list(enumerate(text.split("\n")))

        ta = TransformationApplier()
        for s in scops:
            ta.visit(s, None, None, None)

            l0 = s.block_items[0].coord.line - 1
            l1 = s.block_items[-1].coord.line - 1
            generator = c_generator.CGenerator()
            replacement = generator.visit(s)
            replace_code_between_lines(indexed_lines, replacement, l0, l1)

        # for_node = s.block_items[1].block_items[7]
        # for_node = for_node.stmt.stmt

        lines = [il[1] for il in indexed_lines]
        d.write("\n".join(lines))

if __name__ == "__main__":
    if len(sys.argv) > 2:
        transform(sys.argv[1], sys.argv[2])
    else:
        print("Usage: ./transform.py input-file output-file")

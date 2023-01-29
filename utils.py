import sys
from pycparser import parse_file, c_generator
from pycparser.c_ast import NodeVisitor
import pycparser.c_ast as c_ast
from itertools import zip_longest
from copy import deepcopy
import re

###############################################################################
#                                   generic                                   #
###############################################################################

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

###############################################################################
#                                  find scops                                 #
###############################################################################

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

class IdReplacer(BetterNodeVisitor):
    def __init__(self, rules):
        self.rules = rules

    def visit_ID(self, node, parent, name, index):
        target = self.rules.get(node.name, None)
        if target is not None:
            setattr(parent, name, deepcopy(target))

def dump(node):
    generator = c_generator.CGenerator()
    print(generator.visit(node))

###############################################################################
#                                    unroll                                   #
###############################################################################

# duplicate body and increase step size
def unroll(for_node, n):
    init = for_node.init
    cond = for_node.cond
    iter_id = for_node.next.expr # ASSUMPTION
    body = for_node.stmt

    new_inc = c_ast.BinaryOp("+=", iter_id, c_ast.Constant(type='int', value=n))

    unrolled_body_elements = []
    for i in range(n):
        ir = IdReplacer({iter_id.name: c_ast.BinaryOp("+", iter_id, c_ast.Constant(type="int", value=i))})
        e = deepcopy(body)
        ir.visit(e, None, None, None)
        unrolled_body_elements.append(e)

    new_body = c_ast.Compound(unrolled_body_elements)

    return c_ast.For(init, cond, new_inc, new_body)

def all_equal(l):
    hashed_list = list(map(lambda x: hash(repr(x)), l))
    first_item = hashed_list[0]
    return all([item == first_item for item in hashed_list])

# take a list of for loops, assuming they have the same bounds, make a new for
# loop with their bodies combined
def jam(for_nodes):
    assert all_equal([f.init for f in for_nodes])
    assert all_equal([f.cond for f in for_nodes])
    assert all_equal([f.next for f in for_nodes])
    init = for_nodes[0].init
    cond = for_nodes[0].cond
    inc = for_nodes[0].next
    new_body = c_ast.Compound([f.stmt for f in for_nodes])
    return c_ast.For(init, cond, inc, new_body)

# unroll the top loop and combine the loops that are in the unrolled body
def unroll_and_jam(for_node, n):
    unrolled_for = unroll(for_node, n)
    non_loop_children = list(filter(lambda c: not isinstance(c, c_ast.For), unrolled_for.stmt.block_items))
    loop_children = list(filter(lambda c: isinstance(c, c_ast.For), unrolled_for.stmt.block_items))
    jammed_loop_child = jam(loop_children)
    new_children = non_loop_children + [jammed_loop_child]
    unrolled_for.stmt = c_ast.Compound(new_children)
    return unrolled_for

# traverse the for_node and unroll_and_jam by the amount specified for its
# iterator's name
def unroll_and_jam_loop_nest(for_node, unroll_guide):
    return None

###############################################################################
#                                     licm                                    #
###############################################################################

# checks if node references the for_node iterator at all
def is_loop_invariant(node, for_node):
    return None

# traverses body for invariant statements and expressions, replaces them in the
# loop with id("invariant_i"), and associates that id with the expression in the
# dict
# SIMPLIFICATION: ONLY CONSIDER ARRAY REFERENCES
def replace_invariant_code(for_node, i=0):
    return None, None

# create a compound stmt containing initialization code for loop invariant
# variables and the for loop with the invariant code replaced.
def hoist(for_node, i=0):
    return for_node_comp

# does a post-order traversal of the for loop nest. at each for_loop, hoist
# it. hoisting a higher loop can hoist the initialization created from hoisting
# a lower loop.
def preform_all_hoisting(for_node):
    i = 0
    return None

###############################################################################
#                                 prefetching                                 #
###############################################################################

# traverses body for invariant statements and expressions, replaces them in the
# loop with id("prefetch_j"), and associates that id with the expression in the
# dict.
# SIMPLIFICATION: ONLY CONSIDER ARRAY REFERENCES
def replace_prefetched_data(for_node, i=0):
    return None, None

# create a compound stmt containing initialization code for prefetched data, the
# for loop with the prefetched data replaced, prefetch variable updates at the
# end of the loop, and finally an epilogue which repeats the loop body with
# prefetched data replaced.
def prefetch(for_node, i=0):
    return for_node_comp

# does a post-order traversal of the for loop nest. at each for_loop, add data
# prefetching. there should be no interference between loops in this operation
def preform_all_prefetching(for_node_compound_stmt):
    i = 0
    return None

# ###############################################################################
# #                                  coalescing                                 #
# ###############################################################################

# # from a for loop (i=i0, i<N) which contains an inner for_loop (j=j0, j<M), make
# # a new for loop which traverses their space together (ij = 0, iJ <
# # (N-i0)*(M-j0)). the body immediately recovers the original indices (i = ij /
# # (N-i0) + i0, j = ij % (N-i0) + j0), and then executes each loop's body.
# def coalesce_two_loops(for_node):
#     return None

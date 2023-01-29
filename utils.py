import sys
from pycparser import parse_file, c_generator
from pycparser.c_ast import NodeVisitor
import pycparser.c_ast as c_ast
from itertools import zip_longest
from copy import copy, deepcopy
from functools import reduce
import re

###############################################################################
#                                   generic                                   #
###############################################################################

def grouper(iterable, n, *, incomplete='fill', fillvalue=None):
    "Collect data into non-overlapping fixed-length chunks or blocks"
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

class NodeMapper(object):
    _method_cache = None

    def map(self, node):
        """ Map a node.
        """

        if self._method_cache is None:
            self._method_cache = {}

        mapper = self._method_cache.get(node.__class__.__name__, None)
        if mapper is None:
            method = 'map_' + node.__class__.__name__
            mapper = getattr(self, method, self.generic_map)
            self._method_cache[node.__class__.__name__] = mapper

        return mapper(node)

    def generic_map(self, node):
        """ Called if no explicit visitor function exists for a
            node. Implements preorder visiting of the node.
        """

        if isinstance(node, c_ast.Compound):
            block_items = [self.map(i) for i in node.block_items]
            return c_ast.Compound(block_items)
        else:
            return deepcopy(node)

def loop_lower_bound(for_node):
    if isinstance(for_node.init, c_ast.DeclList):
        return for_node.init.decls[0].init # ASSUMPTION
    elif isinstance(for_node.init, c_ast.Assignment):
        return for_node.init.rvalue

def loop_iter(for_node):
    return for_node.cond.left

def step_amount(for_node):
    inc = for_node.next
    if isinstance(inc, c_ast.UnaryOp):
        if inc.op == "p++":
            return c_ast.Constant(type="int", value=1)
    elif isinstance(inc, c_ast.BinaryOp):
        if inc.op == "+=":
            return inc.right

def concat_nodes(node1, node2):
    if isinstance(node1, c_ast.Compound):
        node1_elements = node1.block_items
        if len(node1_elements) == 0:
            return node2
    elif isinstance(node1, list):
        node1_elements = node1
    else:
        node1_elements = [node1]

    if isinstance(node1, c_ast.Compound):
        node2_elements = node2.block_items
        if len(node1_elements) == 0:
            return node2
    elif isinstance(node2, list):
        node2_elements = node2
    else:
        node2_elements = [node2]

    elements = node1_elements + node2_elements
    return c_ast.Compound(elements)

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
            self.scops.append(c_ast.Compound(scop))

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

    new_body = reduce(concat_nodes, unrolled_body_elements)

    return c_ast.For(init, cond, new_inc, new_body)

def all_equal(l):
    hashed_list = list(map(lambda x: hash(repr(x)), l))
    first_item = hashed_list[0]
    return all([item == first_item for item in hashed_list])

# take a non empty list of for loops, assuming they have the same bounds, make a new for
# loop with their bodies combined
def jam(for_nodes):
    assert len(for_nodes) > 1
    assert all_equal([f.init for f in for_nodes])
    assert all_equal([f.cond for f in for_nodes])
    assert all_equal([f.next for f in for_nodes])
    init = for_nodes[0].init
    cond = for_nodes[0].cond
    inc = for_nodes[0].next
    new_body = reduce(concat_nodes, [f.stmt for f in for_nodes])
    return c_ast.For(init, cond, inc, new_body)

# unroll the top loop and combine the loops that are in the unrolled body
def unroll_and_jam(for_node, n):
    unrolled_for = unroll(for_node, n)
    non_loop_children = list(filter(lambda c: not isinstance(c, c_ast.For), unrolled_for.stmt.block_items))
    loop_children = list(filter(lambda c: isinstance(c, c_ast.For), unrolled_for.stmt.block_items))
    if len(loop_children) > 1:
        jammed_loop_child = jam(loop_children)
        new_children = concat_nodes(non_loop_children, jammed_loop_child)
    else:
        new_children = c_ast.Compound(non_loop_children)

    unrolled_for.stmt = new_children
    return unrolled_for

class Unroller(NodeMapper):
    def __init__(self, unroll_guide):
        self.unroll_guide = unroll_guide

    def map_For(self, node):
        iter_id = loop_iter(node)
        if iter_id.name in self.unroll_guide:
            n = self.unroll_guide[iter_id.name]
            loop = unroll_and_jam(node, n)
        else:
            loop = node

        loop.stmt = self.generic_map(loop.stmt)
        return loop

# traverse the for_node and unroll_and_jam by the amount specified for its
# iterator's name.
def preform_all_unrolling(for_node, unroll_guide):
    u = Unroller(unroll_guide)
    return u.map(for_node)

###############################################################################
#                                     licm                                    #
###############################################################################

class IdFinder(NodeVisitor):
    def __init__(self, target_name):
        self.target_name = target_name
        self.found = False

    def visit_ID(self, node):
        if node.name == self.target_name:
            self.found = True

# checks if node references the for_node iterator at all
def is_loop_invariant(node, iter_id_name):
    id_finder = IdFinder(iter_id_name)
    id_finder.visit(node)
    return not id_finder.found

class InvariantExprReplacer(BetterNodeVisitor):
    def __init__(self, target_name, i):
        self.target_name = target_name
        self.invariant_expr_i = i
        self.name_to_invariant_expr_dict = {}
        self.invariant_expr_to_name_dict = {}

    # Don't recurse into for loop children
    def visit_For(self, node, parent, name, index):
        return

    # Only recurse into right side of assignments
    def visit_Assignment(self, node, parent, name, index):
        self.generic_visit(node.rvalue, parent, name, index)

    def visit_ArrayRef(self, node, parent, name, index):
        if is_loop_invariant(node, self.target_name):
            if hash(repr(node)) not in self.invariant_expr_to_name_dict.keys():
                temp_name = f"invariant_{self.invariant_expr_i}"
                self.invariant_expr_to_name_dict[hash(repr(node))] = temp_name
                self.name_to_invariant_expr_dict[temp_name] = node
                self.invariant_expr_i += 1
            else:
                temp_name = self.invariant_expr_to_name_dict[hash(repr(node))]

            # replace it
            setattr(parent, name, c_ast.ID(name=temp_name))

# traverses body for invariant statements and expressions, replaces them in the
# loop with id("invariant_i"), and associates that id with the expression in a
# dict. returns the for loop with replaced expressions and the dict
# SIMPLIFICATION: ONLY CONSIDER ARRAY REFERENCES
def extract_invariant_code(for_node, i):
    iter_id_name = loop_iter(for_node).name # ASSUMPTION
    ier = InvariantExprReplacer(iter_id_name, i)
    new_body = deepcopy(for_node.stmt)
    ier.visit(new_body, None, None, None)

    init = for_node.init
    cond = for_node.cond
    inc = for_node.next
    return c_ast.For(init, cond, inc, new_body), ier.name_to_invariant_expr_dict

def make_decl(name, init):
    decl = c_ast.Decl(name=name,
                      quals=[],
                      align=[],
                      storage=[],
                      funcspec=[],
                      type=c_ast.TypeDecl(declname=name, quals=[], align=None, type=c_ast.IdentifierType(names=['DATA_TYPE'])),
                      init=init,
                      bitsize=None)

    return decl

# create a compound stmt containing initialization code for loop invariant
# variables and the for loop with the invariant code replaced.
def hoist(for_node, i=0):
    new_for_node, invariants_dict = extract_invariant_code(for_node, i)
    invariant_init = [make_decl(name, init) for name, init in invariants_dict.items()]
    return concat_nodes(invariant_init, new_for_node)

class Hoister(NodeMapper):
    def __init__(self):
        self.invariant_i = 0
        self.applied = False

    def map_For(self, node):
        body = self.map(node.stmt)
        node.stmt = body

        if not self.applied:
            loop = hoist(node, self.invariant_i)
            self.applied = True
        else:
            loop = node

        return loop

# does a post-order traversal of the for loop nest. at each for_loop, hoist
# it. hoisting a higher loop can hoist the initialization created from hoisting
# a lower loop.
def preform_all_hoisting(for_node):
    h = Hoister()
    return h.map(for_node)

###############################################################################
#                                 prefetching                                 #
###############################################################################

class ExprReplacer(BetterNodeVisitor):
    def __init__(self, i):
        self.expr_i = i
        self.name_to_expr_dict = {}
        self.expr_to_name_dict = {}

    # Don't recurse into for loop children
    def visit_For(self, node, parent, name, index):
        return

    # Only recurse into right side of assignments
    def visit_Assignment(self, node, parent, name, index):
        self.generic_visit(node.rvalue, parent, name, index)

    def visit_ArrayRef(self, node, parent, name, index):
        if hash(repr(node)) not in self.expr_to_name_dict.keys():
            temp_name = f"prefetch_{self.expr_i}"
            self.expr_to_name_dict[hash(repr(node))] = temp_name
            self.name_to_expr_dict[temp_name] = node
            self.expr_i += 1
        else:
            temp_name = self.expr_to_name_dict[hash(repr(node))]

        # replace it
        setattr(parent, name, c_ast.ID(name=temp_name))


# traverses body for expressions, replaces them in the
# loop with id("prefetch_j"), and associates that id with the expression in the
# dict.
# SIMPLIFICATION: ONLY CONSIDER ARRAY REFERENCES
def extract_data_for_prefetching(for_node, j):
    iter_id_name = loop_iter(for_node).name # ASSUMPTION
    er = ExprReplacer(j)
    new_body = deepcopy(for_node.stmt)
    er.visit(new_body, None, None, None)

    init = for_node.init
    cond = for_node.cond
    inc = for_node.next
    return c_ast.For(init, cond, inc, new_body), er.name_to_expr_dict

def make_assignment(name, value):
    return c_ast.Assignment(op='=',
                            lvalue=c_ast.ID(name=name),
                            rvalue=value)

# create a compound stmt containing initialization code for prefetched data, the
# for loop with the prefetched data replaced, prefetch variable updates at the
# end of the loop, and finally an epilogue which repeats the loop body with
# prefetched data replaced.
def prefetch(for_node, j=0):
    new_for_node, expr_dict = extract_data_for_prefetching(for_node, j)
    # init = new_for_node.init.rvalue # ASSUMPTION
    iter_id = loop_iter(new_for_node) # ASSUMPTION
    lower_bound = loop_lower_bound(new_for_node)
    inc = step_amount(new_for_node)
    body = new_for_node.stmt

    result_block = []

    # the prologue, which initializes the prefetched expressions
    prefetch_declarations = [deepcopy(make_decl(name, init)) for name, init in expr_dict.items()]
    prologue_block = c_ast.Compound(prefetch_declarations)
    ir = IdReplacer({iter_id.name: lower_bound})
    ir.visit(prologue_block, None, None, None)
    result_block += prologue_block

    # the loop, with prefetching appended to the loop body
    prefetch_assignments = [deepcopy(make_assignment(name, init)) for name, init in expr_dict.items()]
    step_block = c_ast.Compound(prefetch_assignments)
    ir = IdReplacer({iter_id.name: c_ast.BinaryOp("+", iter_id, inc)})
    ir.visit(step_block, None, None, None)
    new_for_node.stmt = concat_nodes(new_for_node.stmt, step_block)
    result_block += [new_for_node]

    # the epilogue, which does the body a final time
    result_block += body

    return c_ast.Compound(result_block)

class Prefetcher(NodeMapper):
    def __init__(self):
        self.prefetch_j = 0
        self.applied = False

    def map_For(self, node):
        body = self.map(node.stmt)
        node.stmt = body

        if not self.applied:
            loop = prefetch(node, self.prefetch_j)
            self.applied = True
        else:
            loop = node

        return loop

# does a post-order traversal of the for loop nest. at each for_loop, hoist
# it. hoisting a higher loop can hoist the initialization created from hoisting
# a lower loop.
def preform_all_prefetching(for_node):
    p = Prefetcher()
    return p.map(for_node)
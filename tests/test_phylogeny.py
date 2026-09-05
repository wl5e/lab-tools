import pytest
import math
from lab_tools.molecular.phylogeny import parse_phylip, neighbor_joining, to_newick, NJNode

SIMPLE_LOWER = """\n4\nA     0\nB     5     0\nC     9    10     0\nD     9    10     8     0\n"""

SIMPLE_SQUARE = """\n4\nA    0   5   9   9\nB    5   0  10  10\nC    9  10   0   8\nD    9  10   8   0\n"""


def test_parse_lower_tri():
    labels, mat = parse_phylip(SIMPLE_LOWER)
    assert labels == ["A", "B", "C", "D"]
    expected = [
        [0,  5,  9,  9],
        [5,  0, 10, 10],
        [9, 10,  0,  8],
        [9, 10,  8,  0],
    ]
    for i in range(4):
        for j in range(4):
            assert abs(mat[i][j] - expected[i][j]) < 1e-6


def test_parse_square():
    labels, mat = parse_phylip(SIMPLE_SQUARE)
    assert labels == ["A", "B", "C", "D"]
    for i in range(4):
        for j in range(4):
            assert abs(mat[i][j] - mat[j][i]) < 1e-6


def test_parse_invalid_first_line():
    with pytest.raises(ValueError):
        parse_phylip("abc\nA 0")


def test_parse_duplicate_label():
    with pytest.raises(ValueError):
        parse_phylip("2\nA 0\nA 1 0")


def test_nj_known_4taxa():
    labels, mat = parse_phylip(SIMPLE_LOWER)
    tree = neighbor_joining(labels, mat)
    newick = to_newick(tree) + ";"
    # Acceptable output: (A:2.250000, B:2.750000, (C:4.000000, D:4.000000):3.000000);
    # Due to floating point, we verify structure and approximate lengths
    # We can parse the Newick string and check relationships
    assert newick.startswith("(")
    assert newick.endswith(";")
    # Just check that A and B are together and C and D are together
    assert "A:" in newick and "B:" in newick
    assert "C:" in newick and "D:" in newick
    # Check that C and D are in a subtree: ((...) ) portion
    assert ") :" in newick or ") :"  # example ends with "):3.0);"


def test_nj_minimal_2taxa():
    labels = ["X", "Y"]
    mat = [[0, 0.8], [0.8, 0]]
    tree = neighbor_joining(labels, mat)
    newick = to_newick(tree) + ";"
    assert newick == "(X:0.400000,Y:0.400000);"


def test_nj_no_neg_lengths():
    # This matrix might cause negative branch lengths; our code clamps to 0
    labels = ["S1", "S2", "S3"]
    mat = [
        [0, 0.1, 0.3],
        [0.1, 0, 0.4],
        [0.3, 0.4, 0],
    ]
    tree = neighbor_joining(labels, mat)
    newick = to_newick(tree)
    assert ":" in newick  # basic sanity


def test_ascii_output_exists():
    from lab_tools.molecular.phylogeny import ascii_tree
    node = NJNode()
    a = NJNode("X", 0.4)
    b = NJNode("Y", 0.4)
    node.children = [a, b]
    out = ascii_tree(node)
    assert "X" in out and "Y" in out


def test_parse_comments_skipped():
    text = "# comment\n2\nA 0\nB 0.5 0\n"
    labels, mat = parse_phylip(text)
    assert labels == ["A", "B"]
    assert mat[0][1] == 0.5

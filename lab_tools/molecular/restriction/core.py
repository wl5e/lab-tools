"""Core logic for restriction enzyme digestion."""

import re
from typing import List, Tuple, Dict
from lab_tools.molecular.restriction.enzymes import ENZYMES


def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(comp[c] for c in reversed(seq))


def find_cut_sites(seq: str, enzymes: List[str]) -> List[Tuple[int, str, int, int, str]]:
    """
    Find all cut positions on top strand.
    Returns list of (cut_pos, enzyme, top_cut, bottom_cut, orientation).
    """
    seq = seq.upper()
    if not re.fullmatch(r'^[ACGT]+$', seq):
        raise ValueError("Sequence must contain only A, C, G, T nucleotides.")
    cuts = []
    for name in enzymes:
        if name not in ENZYMES:
            raise ValueError(f"Unknown enzyme: {name}")
        site, tc, bc = ENZYMES[name]
        L = len(site)
        # forward orientation
        for match in re.finditer(site, seq):
            start = match.start()
            cut_top = start + tc
            cuts.append((cut_top, name, tc, bc, '+'))
        rc_site = reverse_complement(site)
        if rc_site != site:  # avoid double counting palindromic sites
            for match in re.finditer(rc_site, seq):
                start = match.start()
                cut_top = start + L - bc
                cuts.append((cut_top, name, tc, bc, '-'))
    cuts.sort(key=lambda x: x[0])
    return cuts


def digest_fragments(sequence: str, enzymes: List[str], circular: bool = False) -> List[Dict]:
    """
    Simulate digestion and return fragment information.
    """
    seq_len = len(sequence)
    cuts = find_cut_sites(sequence, enzymes)
    # Build mapping cut_pos -> overhang type
    overhang_map = {}
    for cut_pos, _, tc, bc, _ in cuts:
        if tc == bc:
            overhang = 'blunt'
        elif tc < bc:
            overhang = "5' overhang"
        else:
            overhang = "3' overhang"
        # in case of multiple cuts at same position (unlikely), first wins
        if cut_pos not in overhang_map:
            overhang_map[cut_pos] = overhang
    # Unique sorted cut positions
    cut_positions = sorted(set(cp for cp, _, _, _, _ in cuts))

    if not cut_positions:
        # no cuts
        return [{
            'start': 1,
            'end': seq_len,
            'length': seq_len,
            'left_enzyme': None,
            'right_enzyme': None,
            'left_overhang': '',
            'right_overhang': ''
        }]

    fragments = []
    if circular:
        n = len(cut_positions)
        for i in range(n):
            left_pos = cut_positions[i]
            right_pos = cut_positions[(i+1) % n]
            left_overhang = overhang_map.get(left_pos, '')
            right_overhang = overhang_map.get(right_pos, '')
            # fragment length
            if i < n - 1:
                length = right_pos - left_pos
                start = left_pos + 1
                end = right_pos
            else:
                # wrap around
                length = (seq_len - left_pos) + right_pos
                start = left_pos + 1
                end = right_pos
            # Find enzyme names from cuts list at each position
            left_enz = _get_enzyme_at_pos(cuts, left_pos)
            right_enz = _get_enzyme_at_pos(cuts, right_pos)
            fragments.append({
                'start': start,
                'end': end,
                'length': length,
                'left_enzyme': left_enz,
                'right_enzyme': right_enz,
                'left_overhang': left_overhang,
                'right_overhang': right_overhang,
            })
    else:
        # Linear
        # Fragment from start to first cut
        first_pos = cut_positions[0]
        fragments.append({
            'start': 1,
            'end': first_pos,
            'length': first_pos,
            'left_enzyme': None,
            'right_enzyme': _get_enzyme_at_pos(cuts, first_pos),
            'left_overhang': '',
            'right_overhang': overhang_map.get(first_pos, ''),
        })
        # Intermediate fragments
        for i in range(len(cut_positions) - 1):
            lpos = cut_positions[i]
            rpos = cut_positions[i+1]
            fragments.append({
                'start': lpos + 1,
                'end': rpos,
                'length': rpos - lpos,
                'left_enzyme': _get_enzyme_at_pos(cuts, lpos),
                'right_enzyme': _get_enzyme_at_pos(cuts, rpos),
                'left_overhang': overhang_map.get(lpos, ''),
                'right_overhang': overhang_map.get(rpos, ''),
            })
        # Last fragment to end
        last_pos = cut_positions[-1]
        fragments.append({
            'start': last_pos + 1,
            'end': seq_len,
            'length': seq_len - last_pos,
            'left_enzyme': _get_enzyme_at_pos(cuts, last_pos),
            'right_enzyme': None,
            'left_overhang': overhang_map.get(last_pos, ''),
            'right_overhang': '',
        })
    return fragments


def _get_enzyme_at_pos(cuts: List[Tuple[int, str, int, int, str]], pos: int) -> str:
    """Return first enzyme name associated with a given cut position."""
    for cut_pos, name, _, _, _ in cuts:
        if cut_pos == pos:
            return name
    return None

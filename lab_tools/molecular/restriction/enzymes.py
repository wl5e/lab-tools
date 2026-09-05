"""Restriction enzyme database."""

ENZYMES = {
    "EcoRI":  ("GAATTC", 1, 5),
    "BamHI":  ("GGATCC", 1, 5),
    "HindIII": ("AAGCTT", 1, 5),
    "PstI":   ("CTGCAG", 5, 1),
    "SmaI":   ("CCCGGG", 3, 3),
    "XhoI":   ("CTCGAG", 1, 5),
    "NotI":   ("GCGGCCGC", 2, 6),
    "KpnI":   ("GGTACC", 5, 1),
    "SacI":   ("GAGCTC", 5, 1),
    "SpeI":   ("ACTAGT", 1, 5),
}

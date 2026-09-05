import pytest
import math

from lab_tools.molecular.primers import wallace_tm, santalucia_tm, validate_sequence


class TestValidateSequence:
    def test_valid_upper(self):
        assert validate_sequence('ACGTACGT') == 'ACGTACGT'

    def test_valid_lower(self):
        assert validate_sequence('acgt') == 'ACGT'

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="non‑empty"):
            validate_sequence('')

    def test_invalid_characters(self):
        with pytest.raises(ValueError, match="Invalid nucleotides"):
            validate_sequence('ACGU')

    def test_numbers_raise(self):
        with pytest.raises(ValueError):
            validate_sequence('A123')


class TestWallaceTm:
    @pytest.mark.parametrize("seq, expected", [
        ("ATGCATGCATGCAT", 40.0),  # 8 A/T, 6 G/C → 2*8 + 4*6 = 40
        ("AAAA", 8.0),            # 4 A → 2*4 = 8
        ("GGGG", 16.0),
        ("AT", 4.0),
        ("GC", 8.0),
    ])
    def test_wallace_values(self, seq, expected):
        assert wallace_tm(seq) == pytest.approx(expected)

    def test_long_sequence_warns(self, capsys):
        _ = wallace_tm("ACGT" * 4)  # 16 bases
        captured = capsys.readouterr()
        assert "Warning" in captured.err


class TestSantaLuciaTm:
    def test_short_raises(self):
        with pytest.raises(ValueError, match="at least 3 nucleotides"):
            santalucia_tm("AT")

    def test_valid_calculation(self):
        # checks that output is a float within reasonable bounds
        tm = santalucia_tm("CGTAGCTAGCT", primer_conc_nM=500, salt_conc_mM=50)
        assert isinstance(tm, float)
        # Typical Tm for an 11‑mer around 30‑50°C
        assert 20 < tm < 60

    def test_gc_rich_high_tm(self):
        tm_gc = santalucia_tm("GGCCGGCCGG", primer_conc_nM=500, salt_conc_mM=50)
        tm_at = santalucia_tm("ATATATATAT", primer_conc_nM=500, salt_conc_mM=50)
        # GC rich should melt significantly higher than AT rich
        assert tm_gc > tm_at + 20

    def test_concentration_affects_tm(self):
        tm_low_conc = santalucia_tm("ATCGATCGAT", primer_conc_nM=100, salt_conc_mM=50)
        tm_high_conc = santalucia_tm("ATCGATCGAT", primer_conc_nM=1000, salt_conc_mM=50)
        # Higher oligo concentration raises Tm
        assert tm_high_conc > tm_low_conc

    def test_salt_affects_tm(self):
        tm_low_salt = santalucia_tm("ATCGATCGAT", primer_conc_nM=500, salt_conc_mM=10)
        tm_high_salt = santalucia_tm("ATCGATCGAT", primer_conc_nM=500, salt_conc_mM=200)
        assert tm_high_salt > tm_low_salt

    def test_invalid_sequence_raises(self):
        with pytest.raises(ValueError):
            santalucia_tm("ACGTX")

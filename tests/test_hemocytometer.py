import pytest
from lab_tools.lab_ops.hemocytometer import CellCounter


class TestCellCounterInit:
    def test_valid_input(self):
        cc = CellCounter(100, 20, 4, 2.0)
        assert cc.live == 100
        assert cc.dead == 20
        assert cc.squares == 4
        assert cc.dilution == 2.0

    def test_negative_live_raises(self):
        with pytest.raises(ValueError, match="non-negative"):
            CellCounter(-1, 0, 4)

    def test_negative_dead_raises(self):
        with pytest.raises(ValueError, match="non-negative"):
            CellCounter(0, -1, 4)

    def test_zero_squares_raises(self):
        with pytest.raises(ValueError, match="positive"):
            CellCounter(10, 0, 0)

    def test_negative_squares_raises(self):
        with pytest.raises(ValueError, match="positive"):
            CellCounter(10, 0, -2)

    def test_zero_dilution_raises(self):
        with pytest.raises(ValueError, match="positive"):
            CellCounter(10, 0, 4, 0.0)

    def test_negative_dilution_raises(self):
        with pytest.raises(ValueError, match="positive"):
            CellCounter(10, 0, 4, -1.5)


class TestProperties:
    @pytest.fixture
    def standard_counter(self):
        # 120 live, 30 dead, 4 squares, dilution 2 -> total 150
        return CellCounter(120, 30, 4, 2.0)

    def test_total_cells_counted(self, standard_counter):
        assert standard_counter.total_cells_counted == 150

    def test_viability(self, standard_counter):
        assert standard_counter.viability == pytest.approx(80.0)

    def test_viability_all_dead(self):
        cc = CellCounter(0, 10, 4)
        assert cc.viability == 0.0

    def test_viability_no_cells(self):
        cc = CellCounter(0, 0, 4)
        assert cc.viability == 0.0

    def test_average_per_square(self, standard_counter):
        assert standard_counter.average_per_square == 37.5

    def test_cells_per_ml(self, standard_counter):
        # average 37.5 * dilution 2 * factor 1e4 = 750000
        assert standard_counter.cells_per_ml == 750000.0

    def test_viable_cells_per_ml(self, standard_counter):
        # 750000 * (120/150) = 600000
        assert standard_counter.viable_cells_per_ml == 600000.0

    def test_cells_per_ml_zero_counts(self):
        cc = CellCounter(0, 0, 4)
        assert cc.cells_per_ml == 0.0
        assert cc.viable_cells_per_ml == 0.0


class TestTotalViableCells:
    def test_valid_volume(self):
        cc = CellCounter(120, 30, 4, 2.0)  # 600000 viable cells/mL
        assert cc.total_viable_cells(2.0) == 1_200_000.0

    def test_zero_volume(self):
        cc = CellCounter(120, 30, 4, 2.0)
        assert cc.total_viable_cells(0.0) == 0.0

    def test_negative_volume_raises(self):
        cc = CellCounter(120, 30, 4, 2.0)
        with pytest.raises(ValueError, match="non-negative"):
            cc.total_viable_cells(-1.0)

    def test_no_viable_cells(self):
        cc = CellCounter(0, 20, 4)  # viability 0
        assert cc.total_viable_cells(5.0) == 0.0


class TestVolumeForCells:
    def test_desired_cells(self):
        cc = CellCounter(120, 30, 4, 2.0)  # 600000 viable/mL
        assert cc.volume_for_cells(300000) == 0.5

    def test_desired_cells_zero_raises(self):
        cc = CellCounter(10, 0, 4)
        with pytest.raises(ValueError, match="positive"):
            cc.volume_for_cells(0)

    def test_no_viable(self):
        cc = CellCounter(0, 10, 4)
        assert cc.volume_for_cells(100) == float('inf')


class TestDilutionVolumes:
    def test_basic_dilution(self):
        cc = CellCounter(120, 30, 4, 2.0)  # C1 = 600000
        v1, v2 = cc.dilution_volumes(200000, 5.0)
        # C2=200000, V2=5 -> total viable needed = 1e6, V1 = 1e6/600000 = 1.666...
        assert v1 == pytest.approx(5/3)
        assert v2 == pytest.approx(5.0 - 5/3)

    def test_target_conc_exceeds_current(self):
        cc = CellCounter(120, 30, 4, 2.0)  # 600000
        # You cannot increase concentration by dilution
        with pytest.raises(ValueError, match="exceeds final volume"):
            cc.dilution_volumes(1_000_000, 5.0)

    def test_zero_viable(self):
        cc = CellCounter(0, 10, 4)
        with pytest.raises(ValueError, match="zero viable"):
            cc.dilution_volumes(100000, 5.0)

    def test_zero_final_volume_raises(self):
        cc = CellCounter(10, 0, 4)
        with pytest.raises(ValueError, match="positive"):
            cc.dilution_volumes(100000, 0.0)

    def test_negative_target_conc_raises(self):
        cc = CellCounter(10, 0, 4)
        with pytest.raises(ValueError, match="positive"):
            cc.dilution_volumes(-100000, 5.0)

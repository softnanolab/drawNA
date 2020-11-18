from drawNA.oxdna.strand import Strand, generate_helix
from drawNA.oxdna.nucleotide import Nucleotide
import numpy as np
import pandas as pd


def test_Strand():
    strand = Strand(
        [
            Nucleotide(
                "A",
                np.array([1.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "A",
                np.array([2.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "A",
                np.array([3.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "A",
                np.array([4.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "A",
                np.array([5.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "A",
                np.array([6.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
        ]
    )

    assert strand.sequence == "AAAAAA"
    assert len(strand.nucleotides) == 6
    assert len(strand.dataframe.columns.values) == 10

    strand.sequence = "cccaaa"
    strand_1 = strand.copy()
    strand_1.sequence = "tTtgGg"

    assert strand.sequence == "CCCAAA"
    assert strand_1.sequence == "TTTGGG"
    assert len(strand_1.nucleotides) == 6
    assert len(strand_1.dataframe.columns.values) == 10

    print("Basic Strand: \n", strand)
    # print(strand.dataframe)
    strand.copy()
    return


def test_generate_helix_ss():
    ssDNA_helix = generate_helix(n=40, double=False)

    assert type(ssDNA_helix) == list
    assert len(ssDNA_helix) == 1
    assert len(ssDNA_helix[0]) == 40
    assert len(ssDNA_helix[0].nucleotides) == 40
    assert len(ssDNA_helix[0].dataframe.columns.values) == 10

    assert ssDNA_helix[0].index != 0  # required for .top oxDNA file
    assert ssDNA_helix[0].index == 1  # put here for explanatory purposes
    print("ssDNA: \n", ssDNA_helix[0])
    # print(ssDNA_helix[0].dataframe)


def test_generate_helix_ds():
    dsDNA_helix = generate_helix(n=10, double=True, sequence="AGGGACGATG")

    assert type(dsDNA_helix) == list
    assert len(dsDNA_helix) == 2
    assert dsDNA_helix[0].sequence == "AGGGACGATG"
    # second strand should have reverse polarity, complementary pairs: T/A & C/G
    assert dsDNA_helix[1].sequence == "CATCGTCCCT"
    assert len(dsDNA_helix[0]) == len(dsDNA_helix[1])
    

    print("dsDNA: \n", dsDNA_helix[0], "\n", dsDNA_helix[1])


def test_generate_helix_seq():
    ssDNA_with_short_seq = generate_helix(n=10, sequence="AAA")
    strand = ssDNA_with_short_seq[0]

    assert len(strand.nucleotides) == 10
    assert strand.sequence[0:3] == "AAA"

    long_seq = "AGAT" * 5
    ssDNA_with_long_seq = generate_helix(n=10, sequence=long_seq)
    strand = ssDNA_with_long_seq[0]
    assert len(strand.nucleotides) == 10
    assert strand.sequence[0:11] == long_seq[0:10]

def test_strand_transform():
    strand = Strand(
        [
            Nucleotide(
                "A",
                np.array([1.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "A",
                np.array([2.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "A",
                np.array([3.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "A",
                np.array([4.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "A",
                np.array([5.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
            Nucleotide(
                "A",
                np.array([6.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0, 0.0, 1.0]),
            ),
        ]
    )
    strand.translate(np.array([10., 0., 0.]))
    

    strand.rotate([0., 0., 0., 1.])
    strand.rotate([0., 0., 0.])
    strand.rotate(
        np.array([
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., 0.],
        ])
    )

    strand.transform(np.array([
        [1., 0., 0., 0.],
        [0., 1., 0., 0.],
        [0., 0., 1., 0.],
    ]))
    return


if __name__ == "__main__":

    test_Strand()
    test_generate_helix_ss()
    test_generate_helix_ds()
    test_generate_helix_seq()

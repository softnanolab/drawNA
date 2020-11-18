from drawNA.oxdna import Nucleotide
import numpy as np
import warnings


def test_Nucleotide():
    nucleotide = Nucleotide(
        "A",
        np.array([1.0, 0.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 0.0, 1.0]),
    )

    assert len(nucleotide.series) == 10
    assert (nucleotide.pos_back == [0.6, 0.0, 0.0]).all()
    assert (nucleotide.pos_back_rel == [-0.4, 0.0, 0.0]).all()
    assert (nucleotide.pos_base == [1.4, 0.0, 0.0]).all()
    assert (nucleotide.pos_stack == [1.34, 0.0, 0.0]).all()
    assert (nucleotide._a2 == [0.0, 1.0, 0.0]).all()

    print(nucleotide)
    print(nucleotide.series)
    nucleotide.make_5p('A')
    nucleotide.make_across()
    nucleotide.copy()
    return

def test_Nucleotide_across():
    nucleotide = Nucleotide(
        "A",
        np.array([1.0, 0.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 0.0, 1.0]),
    )

    # memory address of each object is stored
    # in the Nucleotide._across attribute
    new = nucleotide.make_across()
    assert id(new._across) == id(nucleotide)
    assert id(nucleotide._across) == id(new)

    # Should correctly raise some warnings and pass the
    # following assertions
    new_2 = nucleotide.make_across()
    assert id(new._across) != id(nucleotide)
    assert id(nucleotide._across) != id(new)
    assert id(new_2._across) == id(nucleotide)
    assert id(nucleotide._across) == id(new_2)

    # check the actual setting of indices works
    nucleotide.index = 0
    new_2.index = 1
    assert nucleotide.across == new_2.index
    assert new_2.across == nucleotide.index   

def test_Nucleotide_transform():
    nucleotide = Nucleotide(
        "A",
        np.array([1.0, 0.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 0.0, 1.0]),
    )

    nucleotide.translate(np.array([10., 0., 0.]))
    assert all(nucleotide.pos_com == [11., 0., 0.])

    nucleotide.rotate([0., 0., 0., 1.])
    nucleotide.rotate([0., 0., 0.])
    nucleotide.rotate(
        np.array([
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., 0.],
        ])
    )

    nucleotide.transform(np.array([
        [1., 0., 0., 0.],
        [0., 1., 0., 0.],
        [0., 0., 1., 0.],
    ]))

    return
    
if __name__ == "__main__":
    test_Nucleotide()
    test_Nucleotide_across()

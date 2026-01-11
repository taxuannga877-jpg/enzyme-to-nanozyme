"""
M9L18
=====

"""

import numpy as np

from stk._internal.topology_graphs.edge import Edge

from .cage import Cage
from .vertices import LinearVertex, NonLinearVertex


class M9L18(Cage):
    """
    Represents a cage topology graph.

    Unoptimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock(
            smiles='[Pd+2]',
            functional_groups=(
                stk.SingleAtom(stk.Pd(0, charge=2))
                for i in range(4)
            ),
            position_matrix=[[0, 0, 0]],
        )

        bb2 = stk.BuildingBlock(
            smiles='C1=NC=CC(C2=CC=CC(C3=CC=NC=C3)=C2)=C1',
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#6]',
                    bonders=(1, ),
                    deleters=(),
                ),
            ],
        )

        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.M9L18(
                building_blocks=(bb1, bb2),
            ),
        )

        moldoc_display_molecule = molecule.Molecule(
            atoms=(
                molecule.Atom(
                    atomic_number=atom.get_atomic_number(),
                    position=position,
                ) for atom, position in zip(
                    cage.get_atoms(),
                    cage.get_position_matrix(),
                )
            ),
            bonds=(
                molecule.Bond(
                    atom1_id=bond.get_atom1().get_id(),
                    atom2_id=bond.get_atom2().get_id(),
                    order=(
                        1
                        if bond.get_order() == 9
                        else bond.get_order()
                    ),
                ) for bond in cage.get_bonds()
            ),
        )

    :class:`.MCHammer` optimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock(
            smiles='[Pd+2]',
            functional_groups=(
                stk.SingleAtom(stk.Pd(0, charge=2))
                for i in range(4)
            ),
            position_matrix=[[0, 0, 0]],
        )

        bb2 = stk.BuildingBlock(
            smiles='C1=NC=CC(C2=CC=CC(C3=CC=NC=C3)=C2)=C1',
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#6]',
                    bonders=(1, ),
                    deleters=(),
                ),
            ],
        )

        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.M9L18(
                building_blocks=(bb1, bb2),
                optimizer=stk.MCHammer(),
            ),
        )

        moldoc_display_molecule = molecule.Molecule(
            atoms=(
                molecule.Atom(
                    atomic_number=atom.get_atomic_number(),
                    position=position,
                ) for atom, position in zip(
                    cage.get_atoms(),
                    cage.get_position_matrix(),
                )
            ),
            bonds=(
                molecule.Bond(
                    atom1_id=bond.get_atom1().get_id(),
                    atom2_id=bond.get_atom2().get_id(),
                    order=(
                        1
                        if bond.get_order() == 9
                        else bond.get_order()
                    ),
                ) for bond in cage.get_bonds()
            ),
        )

    Metal building blocks with four functional groups are
    required for this topology.

    Ligand building blocks with two functional groups are required for
    this topology.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cage-topology-graph-examples`:
    *Multi-Building Block Cage Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 4-functional groups: 0 to 8
        | 2-functional groups: 9 to 17

    See :class:`.Cage` for more details and examples.

    """

    # Use a factor here, because extracted from a crystal structure.
    _factor = 1 / 5
    _vertex_prototypes = (
        NonLinearVertex(
            0,
            np.array([4.7, 4.7, 0.7]) * _factor,
            use_neighbor_placement=False,
        ),
        LinearVertex(
            1,
            np.array([5.3, 4.0, -2.1]) * _factor,
            use_neighbor_placement=False,
        ),
        LinearVertex(
            2,
            np.array([5.4, 0.8, 1.1]) * _factor,
            use_neighbor_placement=False,
        ),
        LinearVertex(
            3,
            np.array([1.8, 4.2, 3.2]) * _factor,
            use_neighbor_placement=False,
        ),
        LinearVertex(
            4,
            np.array([2.8, 6.4, -0.7]) * _factor,
            use_neighbor_placement=False,
        ),
        NonLinearVertex(
            5,
            np.array([-1.8, 3.1, 4.3]) * _factor,
            use_neighbor_placement=False,
        ),
        LinearVertex(
            6,
            np.array([-4.5, 0.5, 3.2]) * _factor,
            use_neighbor_placement=False,
        ),
        LinearVertex(
            7,
            np.array([-1.4, -0.6, 5.3]) * _factor,
            use_neighbor_placement=False,
        ),
        LinearVertex(
            8,
            np.array([-1.3, 5.3, 1.1]) * _factor,
            use_neighbor_placement=False,
        ),
        NonLinearVertex(
            9,
            np.array([-5.8, -2.9, 1.8]) * _factor,
            use_neighbor_placement=False,
        ),
        LinearVertex(
            10,
            np.array([-5.1, -1.3, -1.7]) * _factor,
            use_neighbor_placement=False,
        ),
        LinearVertex(
            11,
            np.array([-4.2, -4.0, 4.0]) * _factor,
            use_neighbor_placement=False,
        ),
        LinearVertex(
            12,
            np.array([-4.6, -5.2, 0.5]) * _factor,
            use_neighbor_placement=False,
        ),
        NonLinearVertex(
            13,
            np.array([4.7, -2.9, 0.5]) * _factor,
            use_neighbor_placement=False,
        ),
        LinearVertex(
            14,
            np.array([2.2, -3.9, 3.3]) * _factor,
            use_neighbor_placement=False,
        ),
        LinearVertex(
            15,
            np.array([4.9, -0.5, -2.5]) * _factor,
            use_neighbor_placement=False,
        ),
        LinearVertex(
            16,
            np.array([1.7, -5.3, -0.4]) * _factor,
            use_neighbor_placement=False,
        ),
        NonLinearVertex(
            17,
            np.array([4.1, 2.8, -4.5]) * _factor,
            use_neighbor_placement=False,
        ),
        LinearVertex(
            18,
            np.array([2.3, 5.1, -4.2]) * _factor,
            use_neighbor_placement=False,
        ),
        LinearVertex(
            19,
            np.array([0.7, 1.0, -5.4]) * _factor,
            use_neighbor_placement=False,
        ),
        NonLinearVertex(
            20,
            np.array([-2.0, -6.4, -0.4]) * _factor,
            use_neighbor_placement=False,
        ),
        LinearVertex(
            21,
            np.array([-1.6, -6.3, 2.5]) * _factor,
            use_neighbor_placement=False,
        ),
        LinearVertex(
            22,
            np.array([-2.5, -3.7, -3.3]) * _factor,
            use_neighbor_placement=False,
        ),
        NonLinearVertex(
            23,
            np.array([-1.4, -4.5, 4.8]) * _factor,
            use_neighbor_placement=False,
        ),
        NonLinearVertex(
            24,
            np.array([0.3, 6.2, -2.3]) * _factor,
            use_neighbor_placement=False,
        ),
        LinearVertex(
            25,
            np.array([-1.9, 3.5, -3.9]) * _factor,
            use_neighbor_placement=False,
        ),
        NonLinearVertex(
            26,
            np.array([-2.9, -0.2, -4.8]) * _factor,
            use_neighbor_placement=False,
        ),
    )
    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_vertex_prototypes[0],
            vertex2=_vertex_prototypes[1],
        ),
        Edge(
            id=1,
            vertex1=_vertex_prototypes[0],
            vertex2=_vertex_prototypes[2],
        ),
        Edge(
            id=2,
            vertex1=_vertex_prototypes[0],
            vertex2=_vertex_prototypes[3],
        ),
        Edge(
            id=3,
            vertex1=_vertex_prototypes[0],
            vertex2=_vertex_prototypes[4],
        ),
        Edge(
            id=4,
            vertex1=_vertex_prototypes[5],
            vertex2=_vertex_prototypes[6],
        ),
        Edge(
            id=5,
            vertex1=_vertex_prototypes[5],
            vertex2=_vertex_prototypes[7],
        ),
        Edge(
            id=6,
            vertex1=_vertex_prototypes[5],
            vertex2=_vertex_prototypes[8],
        ),
        Edge(
            id=7,
            vertex1=_vertex_prototypes[3],
            vertex2=_vertex_prototypes[5],
        ),
        Edge(
            id=8,
            vertex1=_vertex_prototypes[9],
            vertex2=_vertex_prototypes[10],
        ),
        Edge(
            id=9,
            vertex1=_vertex_prototypes[6],
            vertex2=_vertex_prototypes[9],
        ),
        Edge(
            id=10,
            vertex1=_vertex_prototypes[9],
            vertex2=_vertex_prototypes[11],
        ),
        Edge(
            id=11,
            vertex1=_vertex_prototypes[9],
            vertex2=_vertex_prototypes[12],
        ),
        Edge(
            id=12,
            vertex1=_vertex_prototypes[13],
            vertex2=_vertex_prototypes[14],
        ),
        Edge(
            id=13,
            vertex1=_vertex_prototypes[13],
            vertex2=_vertex_prototypes[15],
        ),
        Edge(
            id=14,
            vertex1=_vertex_prototypes[13],
            vertex2=_vertex_prototypes[16],
        ),
        Edge(
            id=15,
            vertex1=_vertex_prototypes[2],
            vertex2=_vertex_prototypes[13],
        ),
        Edge(
            id=16,
            vertex1=_vertex_prototypes[1],
            vertex2=_vertex_prototypes[17],
        ),
        Edge(
            id=17,
            vertex1=_vertex_prototypes[15],
            vertex2=_vertex_prototypes[17],
        ),
        Edge(
            id=18,
            vertex1=_vertex_prototypes[17],
            vertex2=_vertex_prototypes[18],
        ),
        Edge(
            id=19,
            vertex1=_vertex_prototypes[17],
            vertex2=_vertex_prototypes[19],
        ),
        Edge(
            id=20,
            vertex1=_vertex_prototypes[20],
            vertex2=_vertex_prototypes[21],
        ),
        Edge(
            id=21,
            vertex1=_vertex_prototypes[20],
            vertex2=_vertex_prototypes[22],
        ),
        Edge(
            id=22,
            vertex1=_vertex_prototypes[16],
            vertex2=_vertex_prototypes[20],
        ),
        Edge(
            id=23,
            vertex1=_vertex_prototypes[12],
            vertex2=_vertex_prototypes[20],
        ),
        Edge(
            id=24,
            vertex1=_vertex_prototypes[11],
            vertex2=_vertex_prototypes[23],
        ),
        Edge(
            id=25,
            vertex1=_vertex_prototypes[21],
            vertex2=_vertex_prototypes[23],
        ),
        Edge(
            id=26,
            vertex1=_vertex_prototypes[14],
            vertex2=_vertex_prototypes[23],
        ),
        Edge(
            id=27,
            vertex1=_vertex_prototypes[7],
            vertex2=_vertex_prototypes[23],
        ),
        Edge(
            id=28,
            vertex1=_vertex_prototypes[24],
            vertex2=_vertex_prototypes[25],
        ),
        Edge(
            id=29,
            vertex1=_vertex_prototypes[8],
            vertex2=_vertex_prototypes[24],
        ),
        Edge(
            id=30,
            vertex1=_vertex_prototypes[4],
            vertex2=_vertex_prototypes[24],
        ),
        Edge(
            id=31,
            vertex1=_vertex_prototypes[18],
            vertex2=_vertex_prototypes[24],
        ),
        Edge(
            id=32,
            vertex1=_vertex_prototypes[10],
            vertex2=_vertex_prototypes[26],
        ),
        Edge(
            id=33,
            vertex1=_vertex_prototypes[19],
            vertex2=_vertex_prototypes[26],
        ),
        Edge(
            id=34,
            vertex1=_vertex_prototypes[25],
            vertex2=_vertex_prototypes[26],
        ),
        Edge(
            id=35,
            vertex1=_vertex_prototypes[22],
            vertex2=_vertex_prototypes[26],
        ),
    )
    _num_windows = 2
    _num_window_types = 1

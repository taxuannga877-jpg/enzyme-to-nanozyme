import pytest

import stk

from ....case_data import CaseData
from ...building_blocks import get_linker, get_pd_atom


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.M4L8x(
                    building_blocks={
                        get_pd_atom(): range(4),
                        get_linker(): range(4, 12),
                    },
                    reaction_factory=stk.DativeReactionFactory(
                        reaction_factory=stk.GenericReactionFactory(
                            bond_orders={
                                frozenset(
                                    {
                                        stk.GenericFunctionalGroup,
                                        stk.SingleAtom,
                                    }
                                ): 9,
                            },
                        ),
                    ),
                ),
            ),
            smiles=(
                "[H]C1=C([H])C2C3=C([H])C([H])=[N](->[Pd+2]45<-[N]6=C([H])C("
                "[H])=C(C7=C([H])C([H])=C([H])C(C8=C([H])C([H])=[N](->[Pd+2]"
                "(<-[N]9=C([H])C([H])=C(C%10=C([H])C(C%11=C([H])C([H])=[N]->"
                "4C([H])=C%11[H])=C([H])C([H])=C%10[H])C([H])=C9[H])(<-[N]4="
                "C([H])C([H])=C(C9=C([H])C(C%10=C([H])C([H])=[N]->5C([H])=C%"
                "10[H])=C([H])C([H])=C9[H])C([H])=C4[H])<-[N]4=C([H])C([H])="
                "C(C5=C([H])C(C9=C([H])C([H])=[N](->[Pd+2]%10%11<-[N]%12=C(["
                "H])C([H])=C(C%13=C([H])C(C%14=C([H])C([H])=[N](->[Pd+2](<-["
                "N]%15=C([H])C([H])=C(C(=C1[H])C=2[H])C([H])=C%15[H])(<-[N]1"
                "=C([H])C([H])=C(C2=C([H])C([H])=C([H])C(C%15=C([H])C([H])=["
                "N]->%10C([H])=C%15[H])=C2[H])C([H])=C1[H])<-[N]1=C([H])C([H"
                "])=C(C2=C([H])C([H])=C([H])C(C%10=C([H])C([H])=[N]->%11C([H"
                "])=C%10[H])=C2[H])C([H])=C1[H])C([H])=C%14[H])=C([H])C([H])"
                "=C%13[H])C([H])=C%12[H])C([H])=C9[H])=C([H])C([H])=C5[H])C("
                "[H])=C4[H])C([H])=C8[H])=C7[H])C([H])=C6[H])C([H])=C3[H]"
            ),
            name=name,
        ),
    ),
)
def metal_cage_m4l8x(request) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )

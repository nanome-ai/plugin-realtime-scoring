from nanome.api import structure
from marshmallow import Schema, fields


class ScoringOutputSchema(Schema):
    """Output schema that custom scoring functions must return."""

    # Complex index which contains the ligand
    complex_index = fields.Integer(required=False)

    # Values relevant to the entire ligand.
    # Key Value pairs will be printed in the main menu while scoring.
    aggregate_scores = fields.List(
        fields.Dict(keys=fields.Str(), values=fields.Float()),
        required=False)

    # List of tuples, where the first element is the atom index and the second
    # element is the atom's score.
    atom_scores = fields.List(
        fields.Tuple((fields.Int(), fields.Float())),
        required=True)



def extract_residues_from_complex(comp, residue_list, comp_name=None):
    """Copy comp, and remove all residues that are not part of the binding site."""
    new_comp = structure.Complex()
    new_mol = structure.Molecule()
    new_comp.add_molecule(new_mol)
    new_comp.name = comp_name or f'{comp.name}'
    new_comp.index = -1
    new_comp.position = comp.position
    new_comp.rotation = comp.rotation

    binding_site_residue_indices = [r.index for r in residue_list]
    for ch in comp.chains:
        reses_on_chain = [res for res in ch.residues if res.index in binding_site_residue_indices]
        if reses_on_chain:
            new_ch = structure.Chain()
            new_ch.name = ch.name
            new_ch.residues = reses_on_chain
            new_mol.add_chain(new_ch)
    return new_comp

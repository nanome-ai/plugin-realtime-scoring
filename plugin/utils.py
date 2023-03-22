from nanome.api import structure


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

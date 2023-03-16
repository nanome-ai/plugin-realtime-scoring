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

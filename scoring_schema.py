from marshmallow import Schema, fields


class LigandScoringSchema(Schema):
    """Output schema that custom scoring functions must return."""
    complex_index = fields.Integer(required=False)
    aggregate_scores = fields.List(
        fields.Dict(keys=fields.Str(), values=fields.Float()),
        required=False)
    atom_scores = fields.List(
        fields.Tuple((fields.Int(), fields.Float())),
        required=True)

from chemdataextractor import Document
from chemdataextractor.model.model import Compound
from chemdataextractor.model.base import BaseModel, FloatType, StringType, ModelType
from chemdataextractor.model.units import TemperatureModel
from chemdataextractor.parse import I, join


class CeilingTemperature(TemperatureModel):
    """The ceiling temperature of the polymerisation for the monomer"""
    specifier = StringType(parse_expression=(I('Boiling') + I('Point')).add_action(join), required=False)
    # specifier = StringType(parse_expression=(I('Boiling') + I('Point')).add_action(join), required=True)
 
class CeilingTemperature(TemperatureModel):
    """The ceiling temperature of the polymerisation for the monomer"""
    specifier = StringType(parse_expression=(I('Boiling') + I('Point')).add_action(join), required=False)
    # specifier = StringType(parse_expression=(I('Boiling') + I('Point')).add_action(join), required=True)
    compound = ModelType(Compound)   compound = ModelType(Compound)

class Monomer(BaseModel):
    """Monomer of the polymerisation reaction"""
    specifier = StringType(parse_expression=I("monomer"), required=False)
    monomer = ModelType(Compound)

class Monomer(BaseModel):
    """Monomer of the polymerisation reaction"""
    specifier = StringType(parse_expression=I("initiator"), required=False)
    monomer = ModelType(Compound)


class Polymerisation(BaseModel):
    monomer = ModelType(Monomer)
    initiator = ModelType(Initiator)

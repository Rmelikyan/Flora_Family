import datetime
from src.data_collection.ParticleProfile import ParticleProfile
from ctypes import c_uint
def my_converter(o):
    '''
    custumized converter for writing internal data to json
    '''
    if isinstance(o, datetime.datetime):
        return o.__str__()
    if isinstance(o, ParticleProfile):
        return o.__dict__
    if isinstance(o, c_uint):
        return o.value
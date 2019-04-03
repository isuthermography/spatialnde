# Not yet implemented

class parameterization(object):
    
    def __init__(self,**kwargs):
        for kwarg in kwargs:
            if not hasattr(self,kwarg):
                raise AttributeError("Unknown attribute %s" % (kwarg))
            
        setattr(self,kwarg,kwargs[kwarg])
        pass
      
    @classmethod
    def new_from_params(cls,cadrepr,uvparams_params):
        return None
    pass


from . import serialization

class VRMLSerialization(serialization):
    buf=None
    opened=None
    
    def __init__(self,**kwargs):
        self.UseCnt={}  # dictionary of use counters for named elements

        for kwarg in kwargs:
            if hasattr(self,kwarg):
                setattr(self,kwarg,kwargs[kwarg])
                pass
            else:
                raise IndexError(kwarg)
            pass
        pass

    def finish(self):
        self.buf.write("\n  ]\n}\n") # Close wrapping separator opened in tofileorbuffer()
        
        if self.opened:
            self.buf.close()
            pass
        pass
    

    @classmethod
    def tofileorbuffer(cls,buf):
        if isinstance(buf,basestring):
            buf=open(buf,"w")
            opened=True
            pass
        else:
            opened=False
            pass

        buf.write("#VRML V2.0 utf8\n")

        buf.write("Group {\n  children [\n")
        return cls(buf=buf,opened=opened)
    pass


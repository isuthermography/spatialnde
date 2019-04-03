
class serialization(object):
    UseCnt=None  # Dictionary of def names indicating how often they have been generated
    # largely abstract.... specialized by x3dserialization and vrmlserialization
    pass


# VRML/X3D escaping and quoting of MFString components 
def escape_mfstring(mfstring):
    mfstring_esc = "".join([ char if char != '"' and char != '\\' else '\\'+char for char in mfstring ])
    return '\"%s\"' % (mfstring_esc)

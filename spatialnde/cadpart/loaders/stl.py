import sys
import string
import numpy as np

def skipws(stl_text,pos):
    while stl_text[pos] in string.whitespace and pos < len(stl_text):
        pos+=1
        pass
    return pos

def skipnonws(stl_text,pos):
    while stl_text[pos] not in string.whitespace and pos < len(stl_text):
        pos+=1
        pass
    return pos


def loadstl_read_vertex(stl_text,pos):
    if not(stl_text[pos:(pos+6)]=="vertex"):
        raise ValueError("Vertex is %s" % (stl_text[pos:(pos+6)]))
    pos+=6

    pos=skipws(stl_text,pos)

    # read vx
    endpos=skipnonws(stl_text,pos)
    vx=float(stl_text[pos:endpos])
    pos=skipws(stl_text,endpos)

    # read vy
    endpos=skipnonws(stl_text,pos)
    vy=float(stl_text[pos:endpos])
    pos=skipws(stl_text,endpos)
    
    # read vz
    endpos=skipnonws(stl_text,pos)
    vz=float(stl_text[pos:endpos])
    pos=skipws(stl_text,endpos)

    

    return (pos,(vx,vy,vz))

def loadstl_parse_facet(stl_text,pos):
    if not(stl_text[pos:(pos+5)]=="facet"):
        raise ValueError("Facet is %s" % (stl_text[pos:(pos+5)]))
    pos+=5
    pos=skipws(stl_text,pos)
    
    if not(stl_text[pos:(pos+6)]=="normal"):
        raise ValueError("Normal is %s" % (stl_text[pos:(pos+6)]))
    pos+=6

    pos=skipws(stl_text,pos)

    # read ni
    endpos=skipnonws(stl_text,pos)
    ni=float(stl_text[pos:endpos])
    pos=skipws(stl_text,endpos)

    # read nj
    endpos=skipnonws(stl_text,pos)
    nj=float(stl_text[pos:endpos])
    pos=skipws(stl_text,endpos)

    # read nk
    endpos=skipnonws(stl_text,pos)
    nk=float(stl_text[pos:endpos])
    pos=skipws(stl_text,endpos)

    # read outer
    if not(stl_text[pos:(pos+5)]=="outer"):
        raise ValueError("Outer is %s" % (stl_text[pos:(pos+5)]))
    pos+=5
    pos=skipws(stl_text,pos)

    # read loop
    if not(stl_text[pos:(pos+4)]=="loop"):
        raise ValueError("Loop is %s" % (stl_text[pos:(pos+4)]))
    pos+=4
    pos=skipws(stl_text,pos)

    (pos,v1)=loadstl_read_vertex(stl_text,pos)
    (pos,v2)=loadstl_read_vertex(stl_text,pos)
    (pos,v3)=loadstl_read_vertex(stl_text,pos)
    
    # read endloop
    if not(stl_text[pos:(pos+7)]=="endloop"):
        raise ValueError("endloop is %s" % (stl_text[pos:(pos+7)]))
    pos+=7
    pos=skipws(stl_text,pos)

    # read endfacet
    if not(stl_text[pos:(pos+8)]=="endfacet"):
        raise ValueError("endfacet is %s" % (stl_text[pos:(pos+8)]))
    pos+=8
    pos=skipws(stl_text,pos)

    
    return (pos, ((ni,nj,nk),v1,v2,v3))


def loadstl_continue_text(stl_text):
    # stl_text begins after solid
    pos=0
    # Ignore first line: find first cr or newline
    firstnewline=stl_text.find('\n')
    firstcr=stl_text.find('\r')

    pos=min(firstnewline,firstcr)
    if pos==-1:  # one not found
        pos=max(firstnewline,firstcr)
        pass

    pos=skipws(stl_text,pos)
    facets=[]  # consists of normal, vertex1, vertex2, vertex3
    
    while pos < len(stl_text) and stl_text[pos:pos+8] != "endsolid":
        (pos,facet)=loadstl_parse_facet(stl_text,pos)
        facets.append(facet)
        pos=skipws(stl_text,pos)
        pass
    
    return np.array(facets,dtype='d')


def loadstl_continue_binary(fh,firstheader):
    stl_dtype = np.dtype([
        ('Normals', np.float32, (3,)),
        ('Vertex1', np.float32, (3,)),
        ('Vertex2', np.float32, (3,)),
        ('Vertex3', np.float32, (3,)),
        ('attr', '<i2', (1,)),
    ])

    header=firstheader+fh.read(75)
    numtri = np.fromfile(fh,'<u4',1)[0]

    facets=np.fromfile(fh,stl_dtype,numtri)
    
    n = facets['Normals']
    v1 = facets['Vertex1']
    v2 = facets['Vertex2']
    v3 = facets['Vertex3']
    
    facets = np.hstack(((n[:, np.newaxis, :]),(v1[:, np.newaxis, :]), (v2[:, np.newaxis, :]), (v3[:, np.newaxis, :])))
    
    trailerdata=fh.read()
    if len(trailerdata) != 2:
        sys.stderr.write("Warning: loadstl_continue_binary(): Junk at end of STL file\n")
        pass


    return facets


def loadstl(filename):
    # returns a (n_triangles x 4 x 3) array
    # index #1 represents 0 for normal, then 1,2,3 for vertices
    
    fh=open(filename,"rb")
    firstheader=fh.read(5)
    if firstheader=="solid":
        # text format
        facets=loadstl_continue_text(fh.read())
        pass
    else:
        # assume binary format
        facets=loadstl_continue_binary(fh,firstheader)
        pass
    fh.close()
    return facets



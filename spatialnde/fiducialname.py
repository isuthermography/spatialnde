import random as rd

def fiducialname(leading_letter=None):
    """ Generate a name for a fiducial mark. Names will generally have 
        six characters and be upper case. If a leading letter is provided, 
        then that letter will be used for the first character. If the
        leading letter is a vowel then the name will only have five 
        characters rather than six """

    Cons=['B', 'C', 'D', 'F','G' ,'H' ,'J' ,'K' ,'L' ,'M' ,'N' ,'P' ,'Q','R','S','T','V','W','X','Y','Z']
    Vows=['A','E','I','O','U']


    IDarr=[]
    if leading_letter is not None:
        IDarr.append(leading_letter.upper())
        pass

    for i in range(2):
        
        if not(i==0 and leading_letter is not None):
            IDarr.append(rd.sample(Cons,1)[0])
            pass
        if not(i==0 and leading_letter is not None and leading_letter.upper() in Vows):
            IDarr.append(rd.sample(Vows,1)[0])
            pass

        IDarr.append(rd.sample(Cons,1)[0])
        pass

    ID="".join(IDarr)

    return ID

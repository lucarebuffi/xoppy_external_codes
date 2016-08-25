C +++
C       integer         function        iblank
C
C       purpose         Returns the last non-white spot in the string.
C
C       input           a fortran character string.
C       output          Index of the last non-white char in the string.
C                       If there are no empty spaces in the string, then
C                       the lenght is simply returned.
C       hacker          Mumit Khan
C ---
        integer function iblank (str)
        implicit        none
        character*(*)   str
        integer         index, ilen
c
c if the last character in the declared string isn't a white space, simply
c return the length.
c
        index = 1
        ilen = len (str)
        if (str(ilen:ilen).NE.' ') then
            index = ilen + 1
            goto 20
        endif
c
 10     continue
        if (str(index:index).NE.' ') then
            index = index + 1
            goto 10
        endif
 20     continue
        iblank = index - 1
        return
        end

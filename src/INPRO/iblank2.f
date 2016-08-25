C +++
C       integer         function        iblank
C
C       purpose         Returns the last non-white spot in the string, starting by the end.
C
C       input           a fortran character string.
C       output          
C ---
        integer function iblank2 (str)
        implicit        none
        character*(*)   str
        integer         ilen
c
c
        ilen = len (str)

 10     if (str(ilen:ilen).EQ.' ') then
            ilen = ilen - 1
            goto 10
        endif

		continue
        iblank2 = ilen
        return
        end

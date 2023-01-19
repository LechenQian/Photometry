# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 21:53:16 2022

@author: qianl
"""

def generate(numRows: int):
    biglist = [[1]]
    if numRows == 1:
        return biglist
    else:
        for i in range(2,numRows+1):
            temp = []
            for j in range(i):
                if j in [0,i-1] :
                    temp.append(1)
                else:
                    
                    temp.append(biglist[-1][j]+biglist[-1][j-1])
            biglist.append(temp)
    return biglist

numRows = 3
biglist = generate(numRows=numRows)
print(biglist)

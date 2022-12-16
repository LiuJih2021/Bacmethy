# -*- coding: utf-8 -*-
#wangbaonan 
"""
Created on Fri Mar 12 13:44:45 2021

@author: 86183
"""
from __future__ import division
import sys
 
def converse_2_matrix(data):
    print(data)
    data = data.split("\n")
    data.pop()
    print(data)
    #print(data)
    #print("\n")
    data_join=[]
    
    for i in range(len(data)):
        data_join.append(data[i].split(" "))
        print(data[i].split(" "))
               
    for i in range(len(data_join)):
        for j in range(len(data_join[i])):
            #print(data_join[i]) test
            data_join[i][j] = int(data_join[i][j])
    return data_join

def calcute_ppm(conver_col_row_data):
    sum = 0
    A=0
    G=0
    C=0
    T=0
    #创建二维空列表的方法
    ppm=[[None]*4 for i in range(len(conver_col_row_data))]

    for i in range(len(conver_col_row_data)):
         sum = conver_col_row_data[i][0] + conver_col_row_data[i][1] + conver_col_row_data[i][2] + conver_col_row_data[i][3]
         A=conver_col_row_data[i][0]/sum
         G=conver_col_row_data[i][1]/sum
         C=conver_col_row_data[i][2]/sum
         T=conver_col_row_data[i][3]/sum
         for j in range(len(conver_col_row_data[i])):
             ppm[i][0] = A
             ppm[i][1] = G
             ppm[i][2] = C
             ppm[i][3] = T
    return ppm
# def calcute_ppm(data_join):
#     sum_A = 0
#     sum_G = 0
#     sum_C = 0
#     sum_T = 0
#     A=[]
#     G=[]
#     C=[]
#     T=[]
#     Total=[]
#     for i in range(len(data_join)):
#         for j in range(len(data_join[i])):
#             if i==0:
#                 sum_A = sum_A + data_join[i][j]
#             if i==1:
#                 sum_G = sum_G + data_join[i][j]
#             if i==2:
#                 sum_C = sum_C + data_join[i][j]
#             if i==3:
#                 sum_T = sum_T + data_join[i][j]
#     for i in range(len(data_join)):
#         for j in range(len(data_join[i])):
#             if i==0:
#                 A.append((data_join[i][j])/sum_A)
#             if i==1:
#                 G.append((data_join[i][j])/sum_G)
#             if i==2:
#                 C.append((data_join[i][j])/sum_C)
#             if i==3:
#                 T.append((data_join[i][j])/sum_T)
#     Total = zip(A,G,C,T)
#     print(sum_A,sum_G,sum_C,sum_T)
#     #print(A,G,C,T,sum_A,sum_G)
#     return Total
    
    
    #print (sum_A,sum_G) test
    
def calcute_colnum(data):
    space_num = 0
    data = data.split("\n")
    for i in range(len(data[0])):
        if data[0][i] ==' ':
           space_num = space_num + 1 
    col_num = space_num+1
    return col_num
def transdata(data):  
    data = data.replace("A [ ","")
    data = data.replace("G [ ","")
    data = data.replace("C [ ","")
    data = data.replace("T [ ","")
    data = data.replace(" ]","")

    data = data.splitlines()
    data = data[1:]


    for i in range(len(data)):
        data[i]= data[i] +"\n" 
    return data
if __name__ == "__main__":
    matrixfile_path = sys.argv[1]
    with open(matrixfile_path) as f:
        data = f.read()
    data = transdata(data)
    data = ''.join(data)
    
    data_join = converse_2_matrix(data)
    print("\n")    
    col_num = calcute_colnum(data)
    conver_col_row_data = [[data_join[j][i] for j in range(len(data_join))] for i in range(col_num)]
    ppm = calcute_ppm(conver_col_row_data)
    ppm = list(ppm)
    
    #print(data_join)
    #print(conver_col_row_data)
    nsites = data_join[0][0] + data_join[1][0] + data_join[2][0] + data_join[3][0]
    print(ppm)
    with open("ppm.txt", 'w') as file:
        file.write("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\n\nA 0.25 C 0.25 G 0.25 T 0.25\n\nMOTIF MA0080.1 SPI1\nletter-probability matrix: alength= 4 w= %d nsites= %d E= 0" %(len(data_join[0]),nsites))
        for i in range(len(ppm)):
            file.write('\n')
            for  j in range(len(ppm[i])):
                file.write(str(ppm[i][j]))
                file.write(' ')
        file.write('\n')
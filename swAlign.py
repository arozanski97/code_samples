#!/usr/bin/env python3

# Write your code here
import sys

def traceback(matrix1,numrow,numcol,seq1,seq2,trace):
    #SW rule: If value during traceback equals 0 end trace
    if matrix1[numrow][numcol]==[0][0]:
        trace=[]
        return(trace)
    else:
    #Set up values for comparison
        left=matrix1[numrow][numcol-1]
        top=matrix1[numrow-1][numcol]
        diag=matrix1[numrow-1][numcol-1]
        myList=[]
        if seq1[(numrow-1)] == seq2[(numcol-1)]:
            myList.extend([seq1[numrow-1],seq2[numcol-1]])
            addrow=-1
            addcol=-1
        elif(diag>=left):
            myList.extend([seq1[numrow-1],seq2[numcol-1]])
            addrow=-1
            addcol=-1
        elif(top>=left):
            myList.extend([seq1[numrow-1],"-"])
            addrow=-1
            addcol=0
        else:
            addcol=-1
            addrow=0
            myList.extend(["-",seq2[numcol-1]])
        row=numrow+addrow
        col=numcol+addcol
        trace=traceback(matrix1,row,col,seq1,seq2,trace)
        trace.extend(myList)
        return(trace)
def main():
#Get the path to the seq files                                                                                                                                                                              
    seqOne=str(sys.argv[1])
    seqTwo=str(sys.argv[2])
#Open the files that contain the two sequences                                                                                                                                                              
    f1=open(seqOne,"r")
    f2=open(seqTwo,"r")
#Remove the header from the files                                                                                                                                                                           
    f1.readline()
    f2.readline()
#Strip the new line from each seq and save it to a string                                                                                                                                                   
    seq1r=f1.readlines()
    seq1=seq1r[0].rstrip()
    seq2r=f2.readlines()
    seq2=seq2r[0].rstrip()
#Set up variables                                                                                                                                                                                           
    matrix=[]
    mymax=0
    coords=[]
#Create empty matrix of dimensions in accordance to seq1 and seq2 lengths                                                                                                                                   
    for i in range(len(seq1)+1):
        new=[]
        new1=[]
        for j in range(len(seq2)+1):
            new.append(None)
            new1.append(None)
        matrix.append(new)

#Fill the matrix, the first column and row consist of 0's                                                                                                                                                   
    for i in range(len(matrix)):
        matrix[i][0]=0
        for j in range(len(matrix[i])):
            if i==0:
                matrix[i][j]=0
#If not in the the first row or column calculate the value of the cell based upon rules                                                                                                                     
            elif i!=0 and j!=0:
                if seq1[i-1] == seq2[j-1]:
                    Match=1
                else:
                    Match=-1
                left=matrix[i][j-1]-1
                top=matrix[i-1][j]-1
                diagonal=matrix[i-1][j-1]+Match
                if (left>=top) and (left>=diagonal):
                    largest=left
                elif(top>=left) and (top>=diagonal):
                    largest=top
                else:
                    largest=diagonal
#SW rule: if the value is less than 0 set it to zero                                                                                                                                                        
                if largest < 0:
                    largest =0
                matrix[i][j]=largest
#Keep track of the largest value in matrix for backtracking                                                                                                                                                 
                if matrix[i][j] > mymax:
                    mymax=matrix[i][j]
                    coords=[i,j]


#Set up values for backtracking
            
    trace=[]
    string1=''
    string2=''
    connect=''
    alignmentscore=0
#Call traceback starting at coords of highest value 
    trace=traceback(matrix,coords[0],coords[1],seq1,seq2,trace)
#Split traceback list into its separate strings 
    for i in range(len(trace)):
        if i % 2==0:
            string1+=trace[i]
        else:
            string2+=trace[i]
#Calculate Alignment Score
    for i,j in zip(range(len(string1)),range(len(string2))):
        if string1[i]==string2[i]:
            connect+='|'
            alignmentscore+=1
        elif string1[i]=='-' or string2[i]=='-':
            connect+=" "
            alignmentscore-=1
        else:
            connect+='*'
            alignmentscore-=1
        
    print(string1)
    print(connect)
    print(string2)
    print("Alignment score:",alignmentscore)
    f1.close()
    f2.close()

if __name__=="__main__":
    main()

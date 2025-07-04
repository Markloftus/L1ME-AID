import pandas as pd
import pysam
import os
import ast
import numpy as np
import json
import collections
from tqdm import tqdm
from Bio.Seq import Seq
from Bio import Align
from Bio import SeqIO
import more_itertools as mit
import argparse as ap
from functools import reduce


def find_sequence(sequence, subsequence):
    positions = []
    start = 0
    while True:
        start = sequence.find(subsequence, start)
        if start == -1:
            break
        positions.append(start)
        start += 1
    return positions




def orientationFinder(df):
    tailDict={'SINE/Alu':{'REV':'TTTTGAGA','FORWARD':'TCTCAAAA'}, 'LINE/L1':{'REV':'AAGTTTTAGGG','FORWARD':'CCCTAAAACTT'},'Retroposon/SVA':{'REV1':'TTATTGATC','REV2':'TTATTGATA','FORWARD1':'GATCAATAAA','FORWARD2':'TATCAATAAA'}}
    df2 = df.copy()
    df2['Orientation']='NONE'

    for row in df2.index:
        #print(row)
        if df2.at[row,'Element_Hits'] == 'NONE':
            #print("NO HITS")
            continue
        else:

            elementAnnotation = str(df2.at[row,'Element_Annotation'])
            designation = str(df2.at[row,'Element_Designation'])
            hitList = ast.literal_eval(str(df2.at[row,'Element_Hits']))
            
            if designation == 'LINE/L1':
                
                appearanceCount='NONE'
                mycount=-1
                appearanceOrderList=[]
                
                sense=0
                antisense=0

                for hit in hitList:                    
                    
                    
                    appearanceOrderList.append(str(hit.split()[10])+"_"+str(hit.split()[8]))
                    
                    mycount+=1
                    if str(hit.split()[10]) == designation:
                        
                        if str(hit.split()[8])=='+':
                            sense+=1
                        else:
                            antisense+=1
                    else:
                        
                        if str(hit.split()[10]) == 'Simple_repeat':
                            appearanceCount=mycount
                        continue
                        
                        
                        
                if appearanceCount == 'NONE':
                    pass
                else:
                    if appearanceCount==0:
                        antisense+=1
                    else:
                        sense+=1
                        
                        
                        
                if str(tailDict[designation]['REV']) in str(df2.at[row,'Sequence']):
                    antisense+=.5
                elif str(tailDict[designation]['FORWARD']) in str(df2.at[row,'Sequence']):
                    sense+=.5
                else:
                    pass
                
                                
                        
                if sense>0 and antisense ==0:
                    df2.at[row,'Orientation']='+'
                elif sense==0 and antisense>0:
                    df2.at[row,'Orientation']='-'
                elif sense>antisense:
                    df2.at[row,'Orientation']='+'
                        
                elif antisense>sense:
                    df2.at[row,'Orientation']='-'
                        
                elif sense == antisense and antisense !=0:
                    df2.at[row,'Orientation']='TWIN_PRIMING'
                
                else:
                    continue
                    
                

            #If this isnt a LINE (so we dont have to worry about twin priming as much)
            else:
                sense=0
                antisense=0

                for hit in hitList:

                    if str(hit.split()[9]) == elementAnnotation:

                        if str(hit.split()[8])=='+':
                            sense+=1
                        else:
                            antisense+=1

                    else:
                        continue
                        
                if designation in tailDict.keys():  
                    for orientTag in tailDict[designation]:
                        if 'REV' in orientTag and orientTag in str(df2.at[row,'Sequence']):
                            antisense+=.5
                        elif 'FORWARD' in orientTag and orientTag in str(df2.at[row,'Sequence']):
                            sense+=.5
                        
                        else:
                            pass
                else:
                    pass

                if sense>0 and antisense ==0:
                    df2.at[row,'Orientation']='+'
                elif sense==0 and antisense>0:
                    df2.at[row,'Orientation']='-'
                elif sense>antisense:
                    df2.at[row,'Orientation']='+'
                elif antisense>sense:
                    df2.at[row,'Orientation']='-'
                else:
                    continue
                
    return(df2)
    


def tailCounter(df):
    
    df2= df.copy()
    df2['FILTER_RESULTS']='Good_Row'
    df2['Tail_Begins']='No_Tail_Detected'
    df2['Tail_Type']='No_Tail_Type'
    df2['Tail_Length']=0
    df2['Tail_Seed_Hits']=0
    
    for row in df2.index:
        

        sequence = str(df2.at[row,'Sequence']).upper()
        
        tTail = 'TTTTT'
        tTailFlag=0
        trunningTotal = 0
        tTailList=[]
        
        if tTail in sequence:
            taillength=0
            findings = [int(x) for x in find_sequence(sequence, tTail)]
            if len(findings)>0:
                tTailFlag=1
                tmin = min(findings)
                
                iterable = findings
                df2.at[row,'Tail_Seed_Hits']=len(findings)
                groupings = [list(group) for group in mit.consecutive_groups(iterable)]
                flag=0
                for group in groupings:
                    if flag==0:

                        taillength+= len(tTail) + (len(group)-1)
                        groupEnd = max(group)+len(tTail)
                        flag=1

                    else:

                        if min(group)-groupEnd<=4:
                            taillength+= (len(tTail) + (len(group)-1) + abs(min(group)-groupEnd))
                            groupEnd = max(group)+len(tTail)
                        else:
                            continue

                tTailLength = taillength
                
            else:
                pass
                
        
        aTail = 'AAAAA'
        aTailFlag=0
        arunningTotal = 0
        aTailList=[]
        if aTail in sequence:
            taillength=0
            reverseSequence = str(Seq(sequence).reverse_complement())
            findings = [int(x) for x in find_sequence(reverseSequence, tTail)]
            df2.at[row,'Tail_Seed_Hits']=len(findings)
            if len(findings)>0:
                aTailFlag=1
                amin = min(findings)
                iterable = findings
                groupings = [list(group) for group in mit.consecutive_groups(iterable)]
                flag=0
                for group in groupings:
                    if flag==0:

                        taillength+= len(aTail) + (len(group)-1)
                        groupEnd = max(group)+len(aTail)
                        flag=1

                    else:

                        if min(group)-groupEnd<=4:
                            taillength+= (len(aTail) + (len(group)-1) + abs(min(group)-groupEnd))
                            groupEnd = max(group)+len(aTail)
                        else:
                            continue

                aTailLength = taillength
                
            else:
                pass


        if tTailFlag>0 and aTailFlag==0:
            
            df2.at[row,'Tail_Begins']=tmin
            df2.at[row,'Tail_Type']='Possible_T-Tail'
            df2.at[row,'Tail_Length']=tTailLength
            if df2.at[row,'Orientation']=='TWIN_PRIMING':
            	df2.at[row,'Orientation']='-'
            else:
            	pass
            
        elif tTailFlag==0 and aTailFlag>0:
            
            df2.at[row,'Tail_Begins']=amin
            df2.at[row,'Tail_Type']='Possible_A-Tail'
            df2.at[row,'Tail_Length']=aTailLength
            if df2.at[row,'Orientation']=='TWIN_PRIMING':
            	df2.at[row,'Orientation']='+'
            else:
            	pass

        
        
        elif aTailFlag>0 and tTailFlag>0:
            
            orientation = str(df2.at[row,'Orientation'])
            
            if tTailLength<aTailLength and amin<tmin:
                df2.at[row,'Tail_Begins']=amin
                df2.at[row,'Tail_Type']='Possible_A-Tail*_and_Possible_T-Tail'
                df2.at[row,'Tail_Length']=aTailLength

                if df2.at[row,'Orientation']=='TWIN_PRIMING':
            	    df2.at[row,'Orientation']='+'
                else:
            	    pass

                
            elif tTailLength>aTailLength and tmin<amin:
                df2.at[row,'Tail_Begins']=tmin
                df2.at[row,'Tail_Type']='Possible_A-Tail_and_Possible_T-Tail*'
                df2.at[row,'Tail_Length']=tTailLength
                if df2.at[row,'Orientation']=='TWIN_PRIMING':
            	    df2.at[row,'Orientation']='-'
                else:
            	    pass
                
            else:

                if orientation !='None':
                    
                    if tmin>=50 and amin<50:
                        df2.at[row,'Tail_Begins']=amin
                        df2.at[row,'Tail_Type']='Possible_A-Tail*_and_Possible_T-Tail'
                        df2.at[row,'Tail_Length']=aTailLength
                        if df2.at[row,'Orientation']=='TWIN_PRIMING':
            	            df2.at[row,'Orientation']='+'
                        else:
            	            pass

                    elif tmin<50 and amin>=50:
                        df2.at[row,'Tail_Begins']=tmin
                        df2.at[row,'Tail_Type']='Possible_A-Tail_and_Possible_T-Tail*'
                        df2.at[row,'Tail_Length']=tTailLength
                        if df2.at[row,'Orientation']=='TWIN_PRIMING':
            	            df2.at[row,'Orientation']='-'
                        else:
            	            pass
                    
                        
                    elif aTailLength<tTailLength and (orientation == '-' or orientation =='TWIN_PRIMING'):
                        df2.at[row,'Tail_Begins']=tmin
                        df2.at[row,'Tail_Type']='Possible_A-Tail_and_Possible_T-Tail*'
                        df2.at[row,'Tail_Length']=tTailLength
                        if df2.at[row,'Orientation']=='TWIN_PRIMING':
            	            df2.at[row,'Orientation']='-'
                        else:
            	            pass

                    elif aTailLength>tTailLength and (orientation == '+' or orientation =='TWIN_PRIMING'):
                        df2.at[row,'Tail_Begins']=amin
                        df2.at[row,'Tail_Type']='Possible_A-Tail*_and_Possible_T-Tail'
                        df2.at[row,'Tail_Length']=aTailLength
                        if df2.at[row,'Orientation']=='TWIN_PRIMING':
            	            df2.at[row,'Orientation']='+'
                        else:
            	            pass
                        
                    elif orientation == '+':
                        df2.at[row,'Tail_Begins']=amin
                        df2.at[row,'Tail_Type']='Possible_A-Tail*_and_Possible_T-Tail'
                        df2.at[row,'Tail_Length']=aTailLength
                        if df2.at[row,'Orientation']=='TWIN_PRIMING':
            	            df2.at[row,'Orientation']='+'
                        else:
            	            pass
                    
                    
                    elif orientation == '-':
                        df2.at[row,'Tail_Begins']=tmin
                        df2.at[row,'Tail_Type']='Possible_A-Tail_and_Possible_T-Tail*'
                        df2.at[row,'Tail_Length']=tTailLength
                        if df2.at[row,'Orientation']=='TWIN_PRIMING':
            	            df2.at[row,'Orientation']='-'
                        else:
            	            pass
                    
                    else:
                        continue
                else:
                    continue
                    
                
        
        else:
            df2.at[row,'FILTER_RESULTS']=['No_Tail_Detected']
            continue
            
    return(df2)
            
   
def aluLinker(df):
    df2 = df.copy()
    
    for row in df2.index:
        
        threePrimeFlags=0
        
        if df2.at[row,'Element_Designation'] == 'SINE/Alu' and str(df2.at[row,'Tail_Type'])!='No_Tail_Type':
            
            tailDesignation = str(df2.at[row,'Tail_Type'])
            tailStartSite = int(df2.at[row,'Tail_Begins'])
            elements = ast.literal_eval(str(df2.at[row,'Element_Hits']))
            
            if len(elements)==1:
                for element in elements:

                    splitElement = element.split()

                    if abs(int(splitElement[12])-tailStartSite)>=120 and abs(int(splitElement[12])-tailStartSite)<=150:
                        threePrimeFlags+=1
                    else:
                        continue
            else:
                pass
                        
            if threePrimeFlags>0:
                df2.at[row,'FILTER_RESULTS']=['Alu_Linker_Region_Warning']
            else:
                continue
            
            
        else:
            continue
    return(df2) 



def finalQuickCheck(df):
    
    df2 = df.copy()
    
    for row in df2.index:
            
            if str(df2.at[row,'Element_Designation']) =='LINE/L1':
                if str(df2.at[row,'Element_Annotation']) in [youngElement for youngElement in str(args.youngLINE).split(",")]:
                    continue
                else:
                    if str(df2.at[row,'FILTER_RESULTS']) == 'Good_Row':
                        df2.at[row,'FILTER_RESULTS']=['OLDER_LINE_SUBFAMILY']
                    else:
                        newList = ast.literal_eval(str(df2.at[row,'FILTER_RESULTS']))
                        newList.append('OLDER_LINE_SUBFAMILY')
                        df2.at[row,'FILTER_RESULTS']=newList
                    
                    
            elif str(df2.at[row,'Element_Designation']) =='SINE/Alu':
                
                aluList = ast.literal_eval(str(df2.at[row,'Element_Hits']))
                if len(aluList)==1:
                    continue
                else:
                    
                    aluHitList=[]
                    for hit in aluList:
                        splitHit = hit.split()
                        if 'Alu' in str(splitHit[9]):
                            aluHitList.append(splitHit[9])
                        else:
                            continue
                            
                    if len(set(aluHitList))>1:
                        if str(df2.at[row,'FILTER_RESULTS']) == 'Good_Row':
                            df2.at[row,'FILTER_RESULTS']=['MULTI-ALU_Hits']
                        else:
                            newList = ast.literal_eval(str(df2.at[row,'FILTER_RESULTS']))
                            newList.append('MULTI-ALU_Hits')
                            df2.at[row,'FILTER_RESULTS']=newList

                    else:
                        continue

    return(df2)


def finalQuickCheck2(df):
    
    df2 = df.copy()
    
    for row in df2.index:
            
            if str(df2.at[row,'Element_Designation']) =='LINE/L1':
                if int(df2.at[row,'Sequence_Length']) >int(str(args.Lengths).split(",")[1]):
                    if str(df2.at[row,'FILTER_RESULTS']) == 'Good_Row':
                        df2.at[row,'FILTER_RESULTS']=['LINE_>Len']
                    else:
                        newList = ast.literal_eval(str(df2.at[row,'FILTER_RESULTS']))
                        newList.append('LINE_>10kLen')
                        df2.at[row,'FILTER_RESULTS']=newList
                else:
                    pass
                    
                if float(df2.at[row,'Element_Divergence']) >float(str(args.Divergences).split(",")[1]):
                    if str(df2.at[row,'FILTER_RESULTS']) == 'Good_Row':
                        df2.at[row,'FILTER_RESULTS']=['LINE_DIVERGENCE']
                    else:
                        newList = ast.literal_eval(str(df2.at[row,'FILTER_RESULTS']))
                        newList.append('LINE_DIVERGENCE')
                        df2.at[row,'FILTER_RESULTS']=newList
                else:
                    continue
                    
                    
            elif str(df2.at[row,'Element_Designation']) =='SINE/Alu':
                
                if int(df2.at[row,'Sequence_Length']) > int(str(args.Lengths).split(",")[0]):
                    if str(df2.at[row,'FILTER_RESULTS']) == 'Good_Row':
                        df2.at[row,'FILTER_RESULTS']=['ALU_>Len']
                    else:
                        newList = ast.literal_eval(str(df2.at[row,'FILTER_RESULTS']))
                        newList.append('ALU_>500Len')
                        df2.at[row,'FILTER_RESULTS']=newList
                else:
                    pass
                    
                if float(df2.at[row,'Element_Divergence']) > float(str(args.Divergences).split(",")[0]):
                    if str(df2.at[row,'FILTER_RESULTS']) == 'Good_Row':
                        df2.at[row,'FILTER_RESULTS']=['ALU_DIVERGENCE']
                    else:
                        newList = ast.literal_eval(str(df2.at[row,'FILTER_RESULTS']))
                        newList.append('ALU_DIVERGENCE')
                        df2.at[row,'FILTER_RESULTS']=newList
                else:
                    continue
                        
                        
            elif str(df2.at[row,'Element_Designation']) =='Retroposon/SVA':
                if int(df2.at[row,'Sequence_Length']) >int(str(args.Lengths).split(",")[2]):
                    if str(df2.at[row,'FILTER_RESULTS']) == 'Good_Row':
                        df2.at[row,'FILTER_RESULTS']=['SVA_>Len']
                    else:
                        newList = ast.literal_eval(str(df2.at[row,'FILTER_RESULTS']))
                        newList.append('SVA_>10kLen')
                        df2.at[row,'FILTER_RESULTS']=newList
                else:
                    pass
                    
                if float(df2.at[row,'Element_Divergence']) >float(str(args.Divergences).split(",")[2]):
                    if str(df2.at[row,'FILTER_RESULTS']) == 'Good_Row':
                        df2.at[row,'FILTER_RESULTS']=['SVA_DIVERGENCE']
                    else:
                        newList = ast.literal_eval(str(df2.at[row,'FILTER_RESULTS']))
                        newList.append('SVA_DIVERGENCE')
                        df2.at[row,'FILTER_RESULTS']=newList
                else:
                    continue
            
            else:
                continue

    return(df2)


def tailCounterCheck(df):
    df2 = df.copy()
    
    for row in df2.index:
        if df2.at[row,'Tail_Type']!='No_Tail_Type':
            flag=0

            if int(df2.at[row,'Tail_Begins'])/int(df2.at[row,'Sequence_Length'])>.5 or int(df2.at[row,'Tail_Begins'])>int(args.maxTailStart):
                flag+=1
            else:
                pass

            if flag ==0:
                continue
            else:
                flagg='Bad_Tail_Position'
                if df2.at[row,'FILTER_RESULTS']!='Good_Row':
                    tempFlagList = ast.literal_eval(str(df2.at[row,'FILTER_RESULTS']))
                    tempFlagList.append(flagg)
                    df2.at[row,'FILTER_RESULTS'] = tempFlagList
                else:
                    df2.at[row,'FILTER_RESULTS'] = [flagg]
        else:
            continue
    
    return(df2)



def cleanDataTEPercentage(df):
    
    df2 = df.copy() 
    df2['Element_Annotation']='No_Element_Annotation'
    df2['Element_Divergence']=0.0
    df2['Element_Proportion']='NONE'
    df2['Element_Percentage']=0.0
    
    #For each locus
    for row in tqdm(df2.index): 
        if df2.at[row,'Element_Hits'] == 'NONE':
            continue

        elif int(df2.at[row,'Sequence_Length'])>50000:
            df2.at[row,'Element_Annotation']='NOT_TESTED_BEYOND_UPPER_SEQ_LENGTH'
            df2.at[row,'Element_Divergence']=100
            continue

        else:


            columnName = 'Sequence'

            tempDict={int(x):0 for x in range(1,len(df2.at[row,columnName])+1)}    
            tempDict2={}
            tempDict3={}
            tempDivDict={}
            umbrellaDict={}

            #All the repeatmasker hits
            teHitList = ast.literal_eval(str(df2.at[row,'Element_Hits']))

            #For hit in repeatmasker list
            for hit in teHitList:

                #Split the hit to get each repeatmasker piece
                splitHit = hit.split()
                focusElement = str(splitHit[9])
                tempDict3[focusElement]=str(splitHit[10])
                
                if focusElement in umbrellaDict:
                    umbrellaDict[focusElement].append((float(splitHit[1])*(int(splitHit[6])-int(splitHit[5]))))
                else:
                    umbrellaDict[focusElement]=[(float(splitHit[1])*(int(splitHit[6])-int(splitHit[5])))]

                

                if focusElement in tempDict2.keys():
                    tempDivDict[focusElement].append(float(splitHit[1]))
                    
                else:
                    tempDivDict[focusElement]=[float(splitHit[1])]
                    tempDict2[focusElement]={x:0 for x in range(1,len(df2.at[row,columnName])+1)}

                for coordinate in range(int(splitHit[5]), int(splitHit[6])+1):
                    tempDict[coordinate]+=1
                    tempDict2[focusElement][coordinate]+=1

            #Percentage of the sequence annotated by repeatmasker as something
            tePercentage = len([x for x in tempDict.values() if x >0])/len(tempDict)
            df2.at[row,'Element_Percentage']=float(tePercentage)

            
            if len(tempDict2)==1:
                mykey = [x for x in tempDict2.keys()][0]
                #which umbrella term
                df2.at[row,'Element_Designation']=str(tempDict3[mykey])
                #what is the sequence proportion made up of this element
                df2.at[row,'Element_Proportion']={mykey:tePercentage}
                df2.at[row,'Element_Annotation']=mykey
                df2.at[row,'Element_Divergence']=np.median(tempDivDict[mykey][0])

            else:

                #Figure out the proportion made by this element
                tempDict4 = {x:len([y for y in tempDict2[x].values() if y>0])/len(tempDict) for x in tempDict2.keys()}
                maxKey = max(tempDict4, key=tempDict4.get)
                df2.at[row,'Element_Designation']=tempDict3[str(maxKey)]
                df2.at[row,'Element_Proportion']=tempDict4
                df2.at[row,'Element_Annotation']=maxKey
                df2.at[row,'Element_Divergence']=sum(umbrellaDict[maxKey])/len(tempDict)
        
    return(df2)



def repeatmaskerPatternFilter(df):
    df2 = df.copy()
    annotationList=[]
    df2['Unique_Element_Count']='One_Element'
    for row in df2.index:
        if df2.at[row,'Element_Proportion'] == 'NONE':
            continue
        else:
            temporaryHitList = ast.literal_eval(str(df2.at[row,'Element_Hits']))

            #I dont really care about simple repeats being found with MEIs because they could be tails or TSDs so ignore them
            if df2.at[row,'Element_Designation'] !='Simple_repeat':
                tempDesignation = [x.split()[10] for x in temporaryHitList if x.split()[10] != 'Simple_repeat']
                try:
                    tempNumberings = [x.split()[14] for x in temporaryHitList if x.split()[10] != 'Simple_repeat']
                except:
                    tempNumberings = [x.split()[9] for x in temporaryHitList if x.split()[10] != 'Simple_repeat']

            else:
                tempDesignation = [x.split()[10] for x in temporaryHitList]
                try:
                    tempNumberings = [x.split()[14] for x in temporaryHitList]
                except:
                    tempNumberings=[x.split()[9] for x in temporaryHitList]
                
                
            if len(set(tempDesignation))==1 and len(set(tempNumberings))==1:
                continue
                
            elif len(set(tempDesignation))==1 and len(set(tempNumberings))>1:
                df2.at[row,'Unique_Element_Count'] = 'One_Element_Type_MultiDesignations'
            
            
            elif len(set(tempDesignation))>1 and len(set(tempNumberings))>1:
                df2.at[row,'Unique_Element_Count'] = 'Multiple_Element_Types_MultiDesignations'
            
            
            else:
                print(row,"SOMETHING ODD")
                
    return(df2)


def findTwinPriming(df):
    df2 = df.copy()
    df2['Twin_Priming_Flag']='NONE'
    
    for row in df2.index:
        if df2.at[row,'Element_Designation']=='LINE/L1':
            allhits = ast.literal_eval(str(df2.at[row,'Element_Hits']))
            orientations = []
            for hit in allhits:
                splitHit=hit.split()
                if str(splitHit[10]) == 'LINE/L1':
                    orientations.append(splitHit[8])
                else:
                    continue
            if len(set(orientations))>1:
                df2.at[row,'Twin_Priming_Flag']='FLAG'
            else:
                continue
            
            
        else:
            continue
    return(df2)


#This function checks if the marjority of a sequence was simple repeat but if there is a second annotation that is Alu,L1,SVA then assign that (cases where a tail might make up the majority of the sequence)
def simpleRepeatCheck(df):
    df2 = df.copy()
    
    for row in df2.index:
        if df2.at[row,'Element_Designation']=='Simple_repeat':
            tempDict = {x:float(y) for x,y in ast.literal_eval(str(df2.at[row,'Element_Proportion'])).items() if ')n' not in x}
            if len(tempDict)>0:
                maxKey = str(max(tempDict, key=tempDict.get)).upper()
                if 'ALU' in maxKey and len(tempDict) == 1:
                    #print(row)
                    df2.at[row,'Element_Designation']= 'SINE/Alu'
                    df2.at[row,'Element_Annotation']= str(max(tempDict, key=tempDict.get))

                elif'L1' in maxKey and len(tempDict) == 1:
                    #print(row)
                    df2.at[row,'Element_Designation']= 'LINE/L1'
                    df2.at[row,'Element_Annotation']= str(max(tempDict, key=tempDict.get))

                elif 'SVA' in maxKey and len(tempDict) == 1:
                    #print(row)
                    df2.at[row,'Element_Designation']= 'Retroposon/SVA'
                    df2.at[row,'Element_Annotation']= str(max(tempDict, key=tempDict.get))

                else:
                    continue
            else:
                continue
        else:
            continue
    return(df2)


def find_consecutive_sublists(lst):
    result = []
    sublist = []

    for i, num in enumerate(lst):
        if not sublist or num == sublist[-1] + 1: 
            sublist.append(num)
        else:
            if len(sublist) > 1:
                result.append(sublist)
            sublist = [num]

    if len(sublist) > 1:
        result.append(sublist)

    return result


def pullTSD(df,genome):
    tsdDict={}
    df2 = df.copy()
    
    print("Finding TSDs")
    for row in tqdm(df2.index):
        
        # Pull Nucleotides from Insertion Coordinate
        segment = str(df2.at[row,'CHROM'])+":"+str(int(df2.at[row,'POS']))+"-"+str(int(df2.at[row,'POS'])+40)
        coordinateSeq = str(''.join(pysam.faidx(genome, segment).split()[1:])).upper()
        
        #Align the coordinates to the sequence provided
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        alignments = aligner.align(str(df2.at[row,'Sequence'].upper()), coordinateSeq)
        alignment = alignments[0]

        #Build a dataframe
        seq1=[str(x).upper() for x in alignment[0,:]]
        seq2 =[str(y).upper() for y in alignment[1,:]]
        tempDF = pd.DataFrame(data=[seq1,seq2]).copy()

        #Drop the columns that dont agree
        goodColumns=[]
        columnInfo={}
        for column in tempDF.columns:
            if len(set(tempDF[column]))==1:
                goodColumns.append(column)
                columnInfo[column]=[x for x in set(tempDF[column])][0]
            else:
                continue
                    
        #Build new dataframe with only shared nucleotides and find the longest stretch
        
        tempDF2 = tempDF[goodColumns].copy()
        sublists ={len(x):x for x in find_consecutive_sublists(tempDF2.columns)}
        if len(sublists)==0:
        	consensus = 'NONE'
        else:
        	longestList = sorted(sublists.keys())[-1]
        	consensus = ''.join([columnInfo[x] for x in sublists[longestList]])
        
        #Build a consensus sequence of the columns nucleotides
        if len(consensus)>=10 and consensus != 'NONE':
            tsdDict[row] = consensus
        else:
            
            # Pull Nucleotides from Insertion Coordinate
            segment = str(df2.at[row,'CHROM'])+":"+str(int(df2.at[row,'POS'])-30)+"-"+str(int(df2.at[row,'POS'])+30)
            coordinateSeq = str(''.join(pysam.faidx(genome, segment).split()[1:])).upper()
            
            mySequenceFront = str(df2.at[row,'Sequence'].upper())[:30]
            mySequenceBack = str(df2.at[row,'Sequence'].upper())[-30:]
            
            kmerSequenceList=[]
            kmerBackList=[]
            kmerFrontList=[]
            
            for i in range(10,40):
                x=0
                while x < (len(coordinateSeq)-(i-1)):
                    kmerSequenceList.append(str(coordinateSeq[x:x+i].upper())) 
                    x+=1
                    
                x=0
                while x < (len(mySequenceFront)-(i-1)):
                    kmerFrontList.append(str(mySequenceFront[x:x+i].upper())) 
                    kmerBackList.append(str(mySequenceBack[x:x+i].upper())) 
                    x+=1
                    
            newList=collections.Counter(list(set(kmerSequenceList))+list(set(kmerBackList))+list(set(kmerFrontList)))
            #print('\n')
            #print(newList)
            newDict = {x:len(x) for x,y in newList.items() if y>1}
            if len(newDict)>0:
                tsdDict[row] = max(newDict, key=newDict.get)
            else:
                tsdDict[row]=consensus

    df2['TSD'] = [tsdDict[x] for x in df2.index]    
    return(df2)


def noTailHailMary(df):
	df2 = df.copy()

	subsetDF = df2[(df2['Tail_Type']=='No_Tail_Type') & (df2['Element_Hits']!='NONE')].copy()

	tailDict={'SINE/Alu':{'REV':'TTTGAGA','FORWARD':'TCTCAAA'},'LINE/L1':{'REV':'AAGTTTTAGGG','FORWARD':'CCCTAAAACTT'},'Retroposon/SVA':{'REV1':'TTATTGATC','REV2':'TTATTGATA','FORWARD1':'GATCAATAA','FORWARD2':'TATCAATAA'}}

	for row in subsetDF.index:

		if subsetDF.at[row,'Element_Designation'] in tailDict.keys():

			rtailFlag=0
			ftailFlag=0
			revmatchedKmerKey='TEMP'
			forwardmatchedKmerKey='TEMP'

			for orientKey in tailDict[subsetDF.at[row,'Element_Designation']]:

				if 'REV' in orientKey:

					if tailDict[subsetDF.at[row,'Element_Designation']][orientKey] in subsetDF.at[row,'Sequence']:
						rtailFlag+=1
						revmatchedKmerKey=tailDict[subsetDF.at[row,'Element_Designation']][orientKey]

					else:
						continue


				else:
					if tailDict[subsetDF.at[row,'Element_Designation']][orientKey] in subsetDF.at[row,'Sequence']:
						ftailFlag+=1
						forwardmatchedKmerKey=tailDict[subsetDF.at[row,'Element_Designation']][orientKey]

					else:
						continue

			if rtailFlag>0 and ftailFlag==0:
				df2.at[row,'Tail_Type'] = 'Possible_MutatedOrTruncated_T-Tail'
				df2.at[row,'Tail_Begins'] = len(df2.at[row,'Sequence'].split(revmatchedKmerKey)[0])
				if df2.at[row,'Orientation']=='TWIN_PRIMING':
					df2.at[row,'Orientation']='-'
				else:
					pass

				if 'No_Tail_Detected' in list(df2.at[row,'FILTER_RESULTS']):

					newListFilter = [x for x in df2.at[row,'FILTER_RESULTS'] if x!='No_Tail_Detected']
					df2.at[row,'FILTER_RESULTS'] = newListFilter

				else:
					pass

			elif rtailFlag==0 and ftailFlag>0:
				df2.at[row,'Tail_Type'] = 'Possible_MutatedOrTruncated_A-Tail'
				df2.at[row,'Tail_Begins'] = len(df2.at[row,'Sequence'].split(forwardmatchedKmerKey)[-1])
				if df2.at[row,'Orientation']=='TWIN_PRIMING':
					df2.at[row,'Orientation']='+'
				else:
					pass

				if 'No_Tail_Detected' in list(df2.at[row,'FILTER_RESULTS']):

					newListFilter = [x for x in df2.at[row,'FILTER_RESULTS'] if x!='No_Tail_Detected']
					df2.at[row,'FILTER_RESULTS'] = newListFilter

				else:
					pass

			elif rtailFlag==0 and ftailFlag==0:
				continue

			else:
				df2.at[row,'Tail_Type'] = 'Possible_MutatedOrTruncated_A-Tail_T-Tail'
				rev_pos = len(df2.at[row,'Sequence'].split(revmatchedKmerKey)[0])
				forw_pos = len(df2.at[row,'Sequence'].split(forwardmatchedKmerKey)[-1])

				# Choose the earlier one
				if rev_pos <= forw_pos:
				    df2.at[row,'Tail_Begins'] = rev_pos
				    tail_type = 'Possible_MutatedOrTruncated_T-Tail'
				else:
				    df2.at[row,'Tail_Begins'] = forw_pos
				    tail_type = 'Possible_MutatedOrTruncated_A-Tail'

				if 'No_Tail_Detected' in list(df2.at[row,'FILTER_RESULTS']):

					newListFilter = [x for x in df2.at[row,'FILTER_RESULTS'] if x!='No_Tail_Detected']
					df2.at[row,'FILTER_RESULTS'] = newListFilter

				else:
					pass

	return(df2)




def main():

	global args 
	# Ask for user inputs
	parser = ap.ArgumentParser()
	parser.add_argument("-i", "--input", dest="input", required=True, help="Please specify the path to the fasta file that you ran RepeatMasker on.")
	parser.add_argument("-r", "--repeats", dest="Repeat", required=True, help="Please specify the path to the RepeatMasker output file.")
	parser.add_argument("-rd", "--RepeatDiv", dest='RepeatDiv', default=20)
	parser.add_argument("-d", "--divergences", dest='Divergences', default='6,15,15')
	parser.add_argument("-l", "--lengths", dest='Lengths', default='500,10000,10000')
	parser.add_argument("-u", "--upper", dest='Upper', default=50000)
	parser.add_argument("-t", "--maxTailStart", dest='maxTailStart', default=50)
	parser.add_argument("-yl", "--youngLINE", dest='youngLINE', default='L1HS,L1P1,L1PA1,L1PA2')
	parser.add_argument("-g", "--referenceGenome", dest='referenceGenome', default='DoNotRun')
	parser.add_argument("-o", "--output", dest='Output', required=True, help="Please specify the output path and filename for the csv file")
	args = parser.parse_args()

	sequences ={}

	fasta_sequences = SeqIO.parse(open(str(args.input)),'fasta')
	for fasta in fasta_sequences:
		name, sequence = fasta.id, str(fasta.seq)
		sequences[name]=sequence
	insDF = pd.DataFrame(data=[[x,y] for x,y in sequences.items()], columns=['ID','Sequence']).set_index("ID")

	rmdirectory = str(args.Repeat)
	callDict={}
	with open(rmdirectory, 'r') as file:
		lines_after_header = file.readlines()[3:]
		for line in lines_after_header:
			goodline = ' '.join(line.split())
			if str(goodline.split(" ")[4]) in callDict.keys() and float(goodline.split(" ")[1])<=float(args.RepeatDiv):
				callDict[str(goodline.split(" ")[4])]['Annotations'].append(goodline)
			elif float(goodline.split(" ")[1])<=float(args.RepeatDiv):
				callDict[str(goodline.split(" ")[4])]={'Annotations':[]}
				callDict[str(goodline.split(" ")[4])]['Annotations'].append(goodline)
			else:
				continue
	file.close()
	insDF['Element_Hits']='NONE'
	for row in insDF.index:
		if row in callDict.keys():
			insDF.at[row,'Element_Hits']=callDict[row]['Annotations']
		else:
			continue

	insDF['Sequence_Length']=[len(insDF.at[x,'Sequence']) for x in insDF.index]


	def process_dataframe(insDF, functions):
		return reduce(lambda df, func: func(df), functions, insDF)

	#processing functions order
	processing_functions = [cleanDataTEPercentage, orientationFinder, tailCounter, aluLinker, tailCounterCheck, finalQuickCheck, finalQuickCheck2, repeatmaskerPatternFilter, findTwinPriming, simpleRepeatCheck, noTailHailMary]

	insDF_Filtered9 = process_dataframe(insDF, processing_functions)

	if str(args.referenceGenome) !='DoNotRun':
		print("You selected to check for TSDs so that is starting (takes around 30sec for 1000 loci). If it fails, check to see if your sequence IDs match what is needed.")
		insDF_Filtered9['CHROM'] = [x.split("-")[0] for x in insDF_Filtered9.index]
		insDF_Filtered9['POS'] = [x.split("-")[1] for x in insDF_Filtered9.index]
		insDF_Filtered10 =  pullTSD(insDF_Filtered9,str(args.referenceGenome))
		print('Writing out the results. Remember to double-check things yourself as this program is in beta!')
		insDF_Filtered10.to_csv(str(args.Output),sep='\t')
	else:
		print('Writing out the results. Remember to double-check things yourself as this program is in beta!')
		insDF_Filtered9.to_csv(str(args.Output),sep='\t')


if __name__ == "__main__":
	main()

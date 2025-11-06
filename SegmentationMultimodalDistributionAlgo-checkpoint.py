from scipy.stats import skew
from scipy.stats import kurtosis
from scipy.stats import jarque_bera
import numpy as np
import pandas as pd
import math

def findIndexOfFirstNaN (mas):
    NanContain=False
    NanFirstIndex = -1
    for i in range(len(mas)):
        if(math.isnan(mas[i])):
            NanContain = True
            NanFirstIndex = i
            break
    return (NanContain, NanFirstIndex)

def findIndexOfFirstNotNaN (mas):
    NotNanContain=False
    NotNanFirstIndex = -1
    for i in range(len(mas)):
        if(math.isnan(mas[i])==False):
            NotNanContain = True
            NotNanFirstIndex = i
            break
    return (NotNanContain,NotNanFirstIndex)

def determineNaN_borders(mas):
    
    OnlyNans = False
    NotNanContain, FirstNotNanIndex = findIndexOfFirstNotNaN(mas)
    if (NotNanContain):
        NanContain, FirstNanIndex = findIndexOfFirstNaN(mas[FirstNotNanIndex:])
    
        if (NanContain):
            borders = (FirstNotNanIndex,FirstNotNanIndex+FirstNanIndex-1)
        else:
            borders = (FirstNotNanIndex,len(mas))
    else:
        OnlyNans = True
        borders = (-1,-1)
    
    return (OnlyNans,borders)


def FoundUniformSubsec(mas, brd=(0,1), epsilon=2):
    borders = brd
    dist= (int)(0.5*(borders[1]-borders[0]))
    borders = (borders[0],borders[0]+dist)
    true_max = true_min = true_av = true_med = true_size = true_kurtosis = true_IQR = true_skew = 0
    def interior(mas, dist=1, epsilon=2):
            
            nonlocal borders
            nonlocal true_max, true_min, true_av, true_med, true_size, true_kurtosis, true_IQR, true_skew
            
            if (dist==0):
                 true_med = true_max = true_min = true_av = mas[borders[0]]
                 true_size = 1
                 true_kurtosis = 0
                 true_skew = 0
                 q75, q25 = np.percentile(mas[borders[0]:borders[1]+1], [75 ,25])
                 true_IQR = q75 - q25
                 borders = (borders[0],borders[0]+1)
                 return
            
            else:
                max_left = max(mas[borders[0]:borders[1]])
                min_left = min(mas[borders[0]:borders[1]])
                if (abs(max_left - min_left)<=epsilon):
                    if (dist>=1):
                        for i in range (dist):
                            max_add = max(mas[borders[0]:borders[1]+dist-i])
                            min_add = min(mas[borders[0]:borders[1]+dist-i])
                            if (abs(max_add - min_add)<=epsilon):
                                borders=(borders[0],borders[1]+dist-i)
                                break
                        true_max = max(mas[borders[0]:borders[1]])
                        true_min = min(mas[borders[0]:borders[1]])
                        true_av = np.mean(mas[borders[0]:borders[1]])
                        true_med = np.median(mas[borders[0]:borders[1]])
                        true_size = borders[1]-borders[0]
                        if (true_av==mas[borders[0]] and true_av==mas[borders[1]-1]):
                            true_kurtosis = 0
                            true_skew = 0
                        else:
                            true_kurtosis = kurtosis(mas[borders[0]:borders[1]])
                            true_skew = skew(mas[borders[0]:borders[1]])
                        
                        q75, q25 = np.percentile(mas[borders[0]:borders[1]], [75 ,25])
                        true_IQR = q75 - q25
                        #print("true borders: ", borders)
                else:
                    borders = (borders[0],borders[0]+dist)
                    dist =(int)(0.5*dist)
                    interior(mas, dist = dist, epsilon = epsilon)
            return    
    interior(mas, dist = dist, epsilon = epsilon)
    return borders, true_max, true_min, true_av, true_med, true_size, true_kurtosis, true_IQR, true_skew

def FoundNotNanSubsec(mas, epsilon=2):
    onlyNan, borders = determineNaN_borders(mas)
    return borders
    
def divideClasses(mas,epsilon = 2):
    NotNanSec_df = pd.DataFrame(data = {'start': [0], 'fin': [0]})
    start_prev = start = 0
    finish = len(mas)
    i=0 
    while (start<=finish):
        #print('!!!!',start,finish,'!!!!')
        brd = FoundNotNanSubsec(mas[start:], epsilon=epsilon)
        NotNanSec_df.loc[i] = [start+brd[0], start +brd[1]]
       # print(start,start +(brd[1]-brd[0]))
        i = i+1
        start = start +brd[1] + 1
    NotNanSec_df['fin'].loc[NotNanSec_df.shape[0]-1]=NotNanSec_df['fin'].loc[NotNanSec_df.shape[0]-1]-1
    
    Result_df = pd.DataFrame(data = {'start': [0], 'fin': [0], 'min': [0], 'max': [0], 'average': [0], 
                                     'median': [0], 'size': [0], 'kurtosis': [0], 'IQR': [0], 'skew': [0]})
    
    j=0
    for i in range (NotNanSec_df.shape[0]):
        start = NotNanSec_df['start'].iloc[i]
        finish = NotNanSec_df['fin'].iloc[i]
        
        while (start<=finish):
            #print("initial start : ", start, "initial finish: ", finish)
            brd, true_max, true_min, true_av, true_med, true_size, true_kurtosis, true_IQR, true_skew = FoundUniformSubsec(mas, brd=(start,finish), epsilon=epsilon)
            #print(brd[0], brd[1], true_min, true_max,true_av, start, finish)
            start = brd[1]
            Result_df.loc[j] = [brd[0], brd[1], true_min, true_max, true_av, true_med, true_size, true_kurtosis, true_IQR, true_skew]
            j=j+1
    
    
    return  Result_df
    
def get_jenks_breaks(data_list, number_class):
    data_list.sort()
    mat1 = []
    for i in range(len(data_list) + 1):
        temp = []
        for j in range(number_class + 1):
            temp.append(0)
        mat1.append(temp)
    mat2 = []
    for i in range(len(data_list) + 1):
        temp = []
        for j in range(number_class + 1):
            temp.append(0)
        mat2.append(temp)
    for i in range(1, number_class + 1):
        mat1[1][i] = 1
        mat2[1][i] = 0
        for j in range(2, len(data_list) + 1):
            mat2[j][i] = float('inf')
    v = 0.0
    for l in range(2, len(data_list) + 1):
        s1 = 0.0
        s2 = 0.0
        w = 0.0
        for m in range(1, l + 1):
            i3 = l - m + 1
            val = float(data_list[i3 - 1])
            s2 += val * val
            s1 += val
            w += 1
            v = s2 - (s1 * s1) / w
            i4 = i3 - 1
            if i4 != 0:
                for j in range(2, number_class + 1):
                    if mat2[l][j] >= (v + mat2[i4][j - 1]):
                        mat1[l][j] = i3
                        mat2[l][j] = v + mat2[i4][j - 1]
        mat1[l][1] = 1
        mat2[l][1] = v
    k = len(data_list)
    kclass = []
    for i in range(number_class + 1):
        kclass.append(min(data_list))
    kclass[number_class] = float(data_list[len(data_list) - 1])
    count_num = number_class
    while count_num >= 2:  # print "rank = " + str(mat1[k][count_num])
        idx = int((mat1[k][count_num]) - 2)
        # print "val = " + str(data_list[idx])
        kclass[count_num - 1] = data_list[idx]
        k = int((mat1[k][count_num] - 1))
        count_num -= 1
    return kclass


    


def fill_classes(df,stat,breaks):
    
    for i in range(1,len(breaks)):
        if (df[stat] <= breaks[i]) & (df[stat] >= breaks[i-1]):
            return i-1
    
    
def segmentationMultimodalDistr (data,bins = 10, classes = 2, stat = 'average'):
  
    h = (np.max(np.abs(data)) - np.min(np.abs(data)))/bins
    T_mas = np.linspace(np.min(np.abs(data))+h,np.max(np.abs(data)),bins-1)
    JB_mas_max = []
    df = pd.DataFrame()
    df['data'] =data
        
    for t in T_mas:
        print(t)
        JB_mas = []
        AdapRes = divideClasses(df['data'], epsilon = t)
        if (stat=='average'):
            av = AdapRes['average'].to_numpy().copy()
        elif (stat=='median'):
            av = AdapRes['median'].to_numpy().copy()
        else: 
            print("Wrong statistical charactistic!!!")
            return -1
        
        print("start_Jenks")
        breaks = get_jenks_breaks(av, classes)
        print("finish_Jenks")
        AdapRes['Jenks_class'] = AdapRes.apply(fill_classes, axis = 1, args =(stat, breaks))
        
        df["Adap_Jenks_Class"] = -100
        
        for i in range(AdapRes.shape[0]):
            df["Adap_Jenks_Class"].iloc[(int)(AdapRes['start'].iloc[i]):(int)(AdapRes['fin'].iloc[i])] = AdapRes['Jenks_class'].iloc[i]

        JB_max = 0
        
        for i in range(classes):
            df_sub = df.loc[df["Adap_Jenks_Class"] == i]
            JB_mas.append(jarque_bera(df_sub['data']).statistic)
        
        JB_mas_max.append(np.max(JB_mas))
    
    jb_min = JB_mas_max[0]
    for jb in JB_mas_max:
        if (np.isnan(jb)==False) and (jb<jb_min):
            jb_min = jb
            
        
    T = T_mas[JB_mas_max.index(jb_min)]
        
    AdapRes = divideClasses(df['data'], epsilon = T)
    av = AdapRes[stat].to_numpy().copy()
        
    breaks = get_jenks_breaks(av, classes)
    AdapRes['Jenks_class'] = AdapRes.apply(fill_classes, axis = 1, args =(stat, breaks))
        
    for i in range(AdapRes.shape[0]):
        df["Adap_Jenks_Class"].iloc[(int)(AdapRes['start'].iloc[i]):(int)(AdapRes['fin'].iloc[i])] = AdapRes['Jenks_class'].iloc[i]
     
    return [T, T_mas, JB_mas_max, df]

def segmentationMultimodalDistrForFixedThresh (data, T, classes = 2, stat = 'average'):
  
    df = pd.DataFrame()
    df['data'] =data
            
    AdapRes = divideClasses(df['data'], epsilon = T)
    av = AdapRes[stat].to_numpy().copy()
        
    breaks = get_jenks_breaks(av, classes)
    AdapRes['Jenks_class'] = AdapRes.apply(fill_classes, axis = 1, args =(stat, breaks))
        
    for i in range(AdapRes.shape[0]):
        df["Adap_Jenks_Class"].iloc[(int)(AdapRes['start'].iloc[i]):(int)(AdapRes['fin'].iloc[i])] = AdapRes['Jenks_class'].iloc[i]
     
    return df
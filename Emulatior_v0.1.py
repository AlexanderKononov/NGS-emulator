import numpy as np
import random
import copy
import sys

N_observations=100
N_subclones=2
exprestion=20
non-abberant_cells=0.2
#N_observ=int(raw_input('Size of emulated data?'))
#N_subclones=int(raw_input('Number of emulated subclones?'))
norm_data=np.zeros((N_observations,4))

change_frame=N_observations/(N_subclones+2)

normal_CNA=np.random.poisson(exprestion, N_observations)
subclon_CNA_data={}
for i in range(N_subclones):
	current_subclone_data=normal_CNA
	current_subclone_data[j]+=i*exprestion for j in range(i*change_frame,(i+1)*change_frame)		
		
	subclon_CNA_data[1]=[]
# DP Emulate
for i in norm_data:
    i[0]=1
    i[3]=1
    candidat=random.randint(1, N_observ*100)
    while candidat in norm_data[:,1]:
        candidat=random.randint(1, N_observ*100)
    i[1]=candidat
    i[2]=int(random.normalvariate(100, 20))
norm_data=norm_data[np.argsort(norm_data[:,1])]

# Subclones CNA Emulate
SC=[]
BSC=[]
VSC=[]
snv_dict={}
change_CN=[0, 0.5, 1.5,2]
for i in range(N_subclones):
    SC.append(copy.deepcopy(norm_data))
    BSC.append(copy.deepcopy(norm_data))
    for j in BSC[i]:
        j[3]=j[2]
        j[2]=0
        a=random.random()
        if a<=0.5:
            j[2]=int(random.normalvariate(j[3]/2, 2))
    coint=0
    evol_case=random.randint(N_observ//10,N_observ//5)
    while coint < evol_case:
        ch_obs=random.randint(0,len(SC[i])-1)
        change=change_CN[random.randint(0,3)]
        x_mean=(SC[i][ch_obs][2]*change)
        SC[i][ch_obs][2]=int(random.normalvariate(x_mean, 2))
        BSC[i][ch_obs][3]=SC[i][ch_obs][2]
        a=random.random()
        if a<=0.5:
            if change==2:
                BSC[i][ch_obs][2]=int(random.normalvariate(BSC[i][ch_obs][3]/2, 2))
            elif change==1.5:
                BSC[i][ch_obs][2]=int(random.normalvariate(BSC[i][ch_obs][3]/3, 2))
            elif change==0:
                BSC[i][ch_obs][2]=0
            elif change==0.5:
                BSC[i][ch_obs][2]=0
                loh_stat=random.random()
                if loh_stat<=0.1:
                    SC[i][ch_obs][3]=int(random.normalvariate(SC[i][ch_obs][3]*2, 2))
                    BSC[i][ch_obs][3]=SC[i][ch_obs][2]
        coint+=1
    VSC.append(BSC[i])
    for j in range(len(VSC[i])):
        mut=random.random()
        if mut <=0.5:
            if VSC[i][j][2]==0:
                VSC[i][j][2]=int(random.normalvariate(VSC[i][j][3]/2,2))
            else:
                VSC[i][j][2]=int(random.normalvariate(VSC[i][j][2]/2,2))
        else:
            VSC[i][j][1]=int(random.normalvariate(VSC[i][j][1],N_observ/10))
            VSC[i][j][2]=int(random.normalvariate(VSC[i][j][2]/2,2))
        
        if VSC[i][j][2]<=0:
            VSC[i][j][2]=1
            
        if VSC[i][j][1] not in snv_dict.keys():
            snv_dict[VSC[i][j][1]]=[[VSC[i][j][2]], VSC[i][j][3]]
        else:
            snv_dict[VSC[i][j][1]][0].append(VSC[i][j][2])
        
print(len(snv_dict))
      
CNA_data=copy.deepcopy(norm_data)
for i in range(len(CNA_data)):
    for j in SC:
        CNA_data[i][2]+=j[i][2]
    CNA_data[i][2]/=(len(SC)+1)
    CNA_data[i][2]=int(CNA_data[i][2])

BAF_data=copy.deepcopy(norm_data)
for i in range(len(BAF_data)):
    for j in BSC:
        BAF_data[i][2]+=j[i][2]
        BAF_data[i][3]+=j[i][3]
    BAF_data[i][2]/=(len(BSC)+1)
    BAF_data[i][2]=int(BAF_data[i][2])
    BAF_data[i][3]/=(len(BSC)+1)
    BAF_data[i][3]=int(BAF_data[i][3])
    
# Output
CNA=open('CNA_Emulate.txt', 'w')
for i in CNA_data:
    CNA.write(str(int(i[0]))+'\t'+str(int(i[1]))+'\t'+str(int(i[2]))+'\t'+str(int(i[3]))+'\n')
CNA.close()

BAF=open('BAF_Emulate.txt', 'w')
for i in BAF_data:
    BAF.write(str(int(i[0]))+'\t'+str(int(i[1]))+'\t'+str(int(i[2]))+'\t'+str(int(i[3]))+'\n')
BAF.close()

SNV=open('SNV_Emulate.txt', 'w')
for i in snv_dict.keys():
    snv_read=int(sum(snv_dict[i][0])/len(snv_dict[i][0]))
    SNV.write('1'+'\t'+str(int(i))+'\t'+str(int(snv_read))+'\t'+str(int(snv_dict[i][1]))+'\n')
SNV.close()


print('-----end-----')
#print(data)


# In[ ]:




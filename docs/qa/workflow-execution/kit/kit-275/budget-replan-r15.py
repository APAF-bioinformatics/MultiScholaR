# KIT-275 G2 at the 1.5 x rule (coordinator's ruling): timeout_minutes = ceil(1.5 x the job's L1-measured compute), compute = data sets x seconds / 4 / 60,
# with the declared continuation. Times: L1's own completed-journey means (run 45756994368). python3 -I budget-replan-r15.py
import math
T={'D3':{'GA':30.6,'GK':29.7,'GSC':34.4,'LA':48.1,'LG':24.7,'LM':27.6,'NO':20.9,'PEG2':19.8,'PEL2':19.8},
   'D4':{'GA':32.9,'GK':32.6,'GSC':36.1,'LA':49.6,'LG':28.4,'LM':30.5,'NO':23.2,'PEG2':28.0,'PEL2':26.4},
   'D6':{'GA':37.3,'GK':36.4,'GSC':38.6,'LA':54.1,'LG':33.0,'LM':33.6,'NO':25.7,'PEG2':31.9,'PEL2':30.0}}
mix={'D':dict(NO=2,LM=3,GK=3,LA=3,GA=3),'K':dict(NO=1,LM=8,LG=3,GK=5),'F':dict(NO=1,LM=1,GK=1,LA=1,GA=1),'S':dict(GSC=3,GK=1),'V':dict(NO=2,LM=3,GK=3,LA=2,GA=2),
     'M':dict(NO=1,LM=3,GK=3),'T':dict(NO=1,LM=1,LG=1,GK=1,LA=1,GA=1),'E':dict(NO=2,PEL2=2,PEG2=2),'ID':dict(NO=1,LM=1,GK=1,LA=1,GA=1)}
def sec(ts,d): return sum(sum(n*T[d][c] for c,n in mix[t].items()) for t in ts)
st={'c':0,'e':0,'big':0}
def j(name,ts,R,d,n,counted=True):
    comp=R*sec(ts,d)/4/60; to=math.ceil(1.5*comp)
    if counted: st['c']+=n*(to+10); st['e']+=n*(comp+13); st['big']=max(st['big'],to+10)
    print('| %s | %s | %d | %d | %.1f | %.1f | %d | %d | %.1f |'%(name,d,n,R,sec(ts,d),comp,to,to+10,comp+13)); return to+10
for d in ('D3','D4','D6'): j('D',['D'],50,d,4)
for d in ('D3','D4','D6'): j('K',['K'],50,d,2)
for d in ('D3','D4','D6'): j('F+S',['F','S'],50,d,2)
j('M',['M'],100,'D6',1); j('ID',['ID'],20,'D4',1)
for d in ('D3','D4','D6'): j('T',['T'],100,d,1)
for d in ('D3','D4','D6'): j('E',['E'],100,d,1)
print('counted jobs (tier V excluded): ceiling',st['c'],'expected',round(st['e'],1),'largest',st['big'])
vj=0
for d in ('D3','D4','D6'): vj+=2*(math.ceil(1.5*50*sec(['V'],d)/4/60)+10)
lv=sum(4*sec(['V'],d) for d in ('D3','D4','D6'))/4/60; lvc=math.ceil(1.5*lv)+10
print('tier V jobs (6, wait for KIT-242): ceiling',vj,'| L1-V ceiling',lvc)
counted=st['c']+vj
sc=40; interim=30; reserve=150; corr=st['big']+40
rest=counted+lvc+sc+interim+reserve+corr
print('remaining job-min:',rest,'= core-h',round(rest*4/60,2),'(counted incl V',counted,', L1-V',lvc,', SC',sc,', interim tier-D read',interim,', reserve',reserve,', correction',corr,')')
used=12.43
print('with used',used,'->',round(rest*4/60+used,2),'; without the interim row',round((rest-interim)*4/60+used,2))

# KIT-275 G2 re-plan on the times measured by smoke send 2 (run 45748396833), completed journeys only. python3 -I budget-replan.py
# Measured: D3 and D6 for NO, LM, LA, GA, GK; D3 for LG, PEL2, PEG2. Scaled: D4 = midpoint of D3 and D6; LG, PEL2, PEG2 at D4 and D6 = their D3 time times (LM at that design / LM at D3);
# GSC = send 1's D3 mean (31.9, completed journeys of that send) scaled the same way. Core-hours = job-minutes x 4 / 60.
import math
t3={'NO':15.7,'LM':19.4,'LG':17.4,'LA':31.7,'GA':20.1,'GK':21.0,'PEL2':16.3,'PEG2':15.9,'GSC':31.9}
t6={'NO':16.6,'LM':21.7,'LA':29.8,'GA':22.5,'GK':22.3}
T={}
for d in ('D3','D4','D6'): T[d]={}
for c,v in t3.items():
    T['D3'][c]=v
    if c in t6:
        T['D6'][c]=t6[c]; T['D4'][c]=(v+t6[c])/2
    else:
        T['D4'][c]=v*((t3['LM']+t6['LM'])/2)/t3['LM']; T['D6'][c]=v*t6['LM']/t3['LM']
meas={'D3':['NO','LM','LA','GA','GK','LG','PEL2','PEG2'],'D6':['NO','LM','LA','GA','GK']}
mix={'D':dict(NO=2,LM=3,GK=3,LA=3,GA=3),'K':dict(NO=1,LM=8,LG=3,GK=5),'F':dict(NO=1,LM=1,GK=1,LA=1,GA=1),'S':dict(GSC=3,GK=1),'V':dict(NO=2,LM=3,GK=3,LA=2,GA=2),
     'M':dict(NO=1,LM=3,GK=3),'T':dict(NO=1,LM=1,LG=1,GK=1,LA=1,GA=1),'E':dict(NO=2,PEL2=2,PEG2=2),'ID':dict(NO=1,LM=1,GK=1,LA=1,GA=1)}
def sec(tier,d): return sum(n*T[d][c] for c,n in mix[tier].items())
def job(tiers,R,d):
    tiers=[tiers] if isinstance(tiers,str) else tiers
    s=sum(sec(t,d) for t in tiers); plan=R*s*1.25/4/60; to=math.ceil(1.5*plan); return s,plan,to,to+10,R*s/4/60+13
rows=[];tot_e=tot_c=0;maxjob=0
def add(name,tiers,R,ds,n):
    global tot_e,tot_c,maxjob
    for d in ds:
        s,p,to,c,e=job(tiers,R,d); tot_e+=e*n; tot_c+=c*n; maxjob=max(maxjob,to+10)
        rows.append((name,d,n,R,round(s,1),round(p,1),to,c,round(e,1)))
add('D','D',50,('D3','D4','D6'),4); add('K','K',50,('D3','D4','D6'),2); add('F+S','FS',50,('D3','D4','D6'),2) if False else None
for d in ('D3','D4','D6'):
    s,p,to,c,e=(lambda ss:(ss,50*ss*1.25/4/60,math.ceil(1.5*50*ss*1.25/4/60),math.ceil(1.5*50*ss*1.25/4/60)+10,50*ss/4/60+13))(sec('F',d)+sec('S',d))
    tot_e+=2*e; tot_c+=2*c; maxjob=max(maxjob,c); rows.append(('F+S',d,2,50,round(sec('F',d)+sec('S',d),1),round(p,1),to,c,round(e,1)))
add('V','V',50,('D3','D4','D6'),2); add('M','M',100,('D6',),1); add('ID','ID',20,('D4',),1); add('T','T',100,('D3','D4','D6'),1); add('E','E',100,('D3','D4','D6'),1)
for r in rows: print(r)
# L1 without V: 3 data sets + E9's second run per design (tiers D,K,F,S,T,E; M at D6; ID at D4)
l1=0
for d in ('D3','D4','D6'):
    ss=sum(sec(t,d) for t in 'DKFSTE')+(sec('M',d) if d=='D6' else 0)+(sec('ID',d) if d=='D4' else 0); l1+=4*ss
l1_plan=l1*1.25/4/60; l1_to=math.ceil(1.5*l1_plan)
lv=sum(4*sec('V',d) for d in ('D3','D4','D6')); lv_plan=lv*1.25/4/60; lv_to=math.ceil(1.5*lv_plan)
print('counted expected min',round(tot_e,1),'ceiling min',tot_c,'largest job ceiling',maxjob)
print('L1',round(l1,1),'s plan',round(l1_plan,1),'timeout',l1_to,'ceiling',l1_to+10,'| L1-V',round(lv,1),round(lv_plan,1),lv_to,lv_to+10)
fixed=30+40+150+(maxjob+40)
total_ceiling=tot_c+(l1_to+10)+(lv_to+10)+fixed
total_expected=tot_e+(l1*0.8/4/60+13)+ (lv/4/60+13)+15+20
print('TOTAL ceiling job-min',round(total_ceiling),'core-h',round(total_ceiling*4/60,2),'| expected job-min',round(total_expected),'core-h',round(total_expected*4/60,2))
print('counted ceiling core-h',round(tot_c*4/60,2))

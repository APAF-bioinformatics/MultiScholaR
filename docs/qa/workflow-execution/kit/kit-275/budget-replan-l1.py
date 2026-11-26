# KIT-275 G2 re-plan on the times measured by smoke send 2 (run 45748396833), completed journeys only. python3 -I budget-replan.py
# Measured: D3 and D6 for NO, LM, LA, GA, GK; D3 for LG, PEL2, PEG2. Scaled: D4 = midpoint of D3 and D6; LG, PEL2, PEG2 at D4 and D6 = their D3 time times (LM at that design / LM at D3);
# GSC = send 1's D3 mean (31.9, completed journeys of that send) scaled the same way. Core-hours = job-minutes x 4 / 60.
import math
# G2: L1's OWN completed-journey means (run 45756994368, 4 workers, mixed tiers), every class at every design, nothing scaled
T={'D3':{'GA':30.6,'GK':29.7,'GSC':34.4,'LA':48.1,'LG':24.7,'LM':27.6,'NO':20.9,'PEG2':19.8,'PEL2':19.8},
   'D4':{'GA':32.9,'GK':32.6,'GSC':36.1,'LA':49.6,'LG':28.4,'LM':30.5,'NO':23.2,'PEG2':28.0,'PEL2':26.4},
   'D6':{'GA':37.3,'GK':36.4,'GSC':38.6,'LA':54.1,'LG':33.0,'LM':33.6,'NO':25.7,'PEG2':31.9,'PEL2':30.0}}
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
total_ceiling=tot_c+(l1_to+25)+(lv_to+10)+fixed   # L1's request: the sim line, the self-test line (15) and 10 of set-up
total_expected=tot_e+(l1/4/60+13)+(lv/4/60+13)+15+20
print('TOTAL ceiling job-min',round(total_ceiling),'core-h',round(total_ceiling*4/60,2),'| expected job-min',round(total_expected),'core-h',round(total_expected*4/60,2))
print('counted ceiling core-h',round(tot_c*4/60,2))

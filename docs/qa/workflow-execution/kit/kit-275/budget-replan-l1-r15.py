# KIT-275 G2 at the 1.5 x rule (coordinator's ruling): timeout_minutes = ceil(1.5 x the job's L1-measured compute), compute = data sets x seconds / 4 / 60,
# together with the declared continuation. Times: L1's own completed-journey means (run 45756994368). python3 -I budget-replan-l1-r15.py
import math
exec(open('budget-replan-l1.py').read().split("def sec(")[0].split("# G2: L1's OWN")[1].join(["# G2\n",""]) if False else "")
T={'D3':{'GA':30.6,'GK':29.7,'GSC':34.4,'LA':48.1,'LG':24.7,'LM':27.6,'NO':20.9,'PEG2':19.8,'PEL2':19.8},
   'D4':{'GA':32.9,'GK':32.6,'GSC':36.1,'LA':49.6,'LG':28.4,'LM':30.5,'NO':23.2,'PEG2':28.0,'PEL2':26.4},
   'D6':{'GA':37.3,'GK':36.4,'GSC':38.6,'LA':54.1,'LG':33.0,'LM':33.6,'NO':25.7,'PEG2':31.9,'PEL2':30.0}}
mix={'D':dict(NO=2,LM=3,GK=3,LA=3,GA=3),'K':dict(NO=1,LM=8,LG=3,GK=5),'F':dict(NO=1,LM=1,GK=1,LA=1,GA=1),'S':dict(GSC=3,GK=1),'V':dict(NO=2,LM=3,GK=3,LA=2,GA=2),
     'M':dict(NO=1,LM=3,GK=3),'T':dict(NO=1,LM=1,LG=1,GK=1,LA=1,GA=1),'E':dict(NO=2,PEL2=2,PEG2=2),'ID':dict(NO=1,LM=1,GK=1,LA=1,GA=1)}
def sec(t,d): return sum(n*T[d][c] for c,n in mix[t].items())
def job(tiers,R,d):
    s=sum(sec(t,d) for t in tiers); comp=R*s/4/60; to=math.ceil(1.5*comp); return s,comp,to,to+10,comp+13
rows=[]; tot_c=tot_e=0; big=0; counted=0
def add(name,tiers,R,ds,n,counted_job=True):
    global tot_c,tot_e,big
    for d in ds:
        s,comp,to,c,e=job(tiers,R,d); tot_c+=c*n; tot_e+=e*n; big=max(big,c)
        rows.append((name,d,n,R,round(s,1),round(comp,1),to,c,round(e,1)))
add('D','D',50,('D3','D4','D6'),4); add('K','K',50,('D3','D4','D6'),2); add('F+S','FS',50,('D3','D4','D6'),2)
add('M','M',100,('D6',),1); add('ID','ID',20,('D4',),1); add('T','T',100,('D3','D4','D6'),1); add('E','E',100,('D3','D4','D6'),1)

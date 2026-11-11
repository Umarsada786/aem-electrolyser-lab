"""Reproduce calibration and fixed chronological holdout evaluation."""
from pathlib import Path
import csv,json
import numpy as np
from scipy.optimize import least_squares
from model import Parameters,predict,F,VM
ROOT=Path(__file__).parent
rows=list(csv.DictReader((ROOT/'data/experimental.csv').open()))
rows=[{k:float(v) for k,v in r.items()} for r in rows]
active=[r for r in rows if r['current_density_a_cm2']>0 and r['voltage_v']>0]
cut=int(len(active)*0.7)
train,test=active[:cut],active[cut:]
def fit(rs):
 j=np.array([r['current_density_a_cm2'] for r in rs]);t=np.array([r['temperature_c'] for r in rs]);p=np.array([r['pressure_bar'] for r in rs]);y=np.array([r['voltage_v'] for r in rs]);tk=t+273.15
 base=1.229-0.0009*(tk-298.15)+8.314462618*tk/(2*F)*np.log(p)
 def fun(x):return base+x[0]*tk/328.15*np.arcsinh(j/(2*np.exp(x[1])))+x[2]*j+x[3]*(t-55)*j-y
 sol=least_squares(fun,[.06,np.log(.001),.3,-.002],bounds=([.001,np.log(1e-8),0,-.03],[.3,np.log(.1),2,.03]),max_nfev=20000)
 area=float(np.median([r['power_kw']*1000/(r['current_density_a_cm2']*r['voltage_v']) for r in rs]))
 eta=float(np.median([r['hydrogen_nm3_h']*2*F/(r['current_density_a_cm2']*area*VM*3600) for r in rs]))
 par=Parameters(sol.x[0],float(np.exp(sol.x[1])),sol.x[2],sol.x[3],area,min(1,eta))
 return par,{'optimizer_success':bool(sol.success),'jacobian_condition_number':float(np.linalg.cond(sol.jac)),'unconstrained_faradaic_efficiency':eta}
par,diagnostics=fit(train)
def prediction(r):return predict(r['current_density_a_cm2'],r['temperature_c'],r['pressure_bar'],par)
def metrics(rs):
 out={}
 for key in ['voltage_v','power_kw','hydrogen_nm3_h']:
  y=np.array([r[key] for r in rs]);pred=np.array([prediction(r)[key] for r in rs]);err=pred-y
  out[key]={'rmse':float(np.sqrt(np.mean(err**2))),'mae':float(np.mean(abs(err))),'bias':float(np.mean(err)),'r2':float(1-sum(err**2)/sum((y-y.mean())**2)) if np.ptp(y)>0 else None,'max_abs_error':float(max(abs(err)))}
 return out
report={'split':'First 70% active rows for calibration; final 30% chronological holdout. No refit on holdout.','total_rows':len(rows),'excluded_zero_current_rows':len(rows)-len(active),'train_n':len(train),'test_n':len(test),'train_last_time_s':train[-1]['time_s'],'test_first_time_s':test[0]['time_s'],'parameters':vars(par),'diagnostics':diagnostics,'calibration':metrics(train),'holdout':metrics(test),'domain':{k:[min(r[k] for r in train),max(r[k] for r in train)] for k in ['current_density_a_cm2','temperature_c','pressure_bar']}}
(ROOT/'docs/validation.json').write_text(json.dumps(report,indent=2))
preds=[]
for r in active:
 z=dict(r);z['split']='calibration' if r['time_s']<=train[-1]['time_s'] else 'holdout';z.update({'pred_'+k:v for k,v in prediction(r).items()});preds.append(z)
with (ROOT/'data/predictions.csv').open('w') as f:
 out=csv.DictWriter(f,fieldnames=preds[0].keys());out.writeheader();out.writerows(preds)
(ROOT/'docs/dataset.js').write_text('window.AEM_DATA='+json.dumps({'report':report,'rows':preds},separators=(',',':'))+';')
print(json.dumps(report,indent=2))

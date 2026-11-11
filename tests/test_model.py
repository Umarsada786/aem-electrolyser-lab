import unittest,math,json,subprocess
from pathlib import Path
from model import predict,Parameters,thermal_step,F
class ModelTests(unittest.TestCase):
 def test_zero_current(self):
  r=predict(0);self.assertEqual(r['hydrogen_kg_h'],0);self.assertEqual(r['power_kw'],0);self.assertIsNone(r['specific_energy_kwh_kg'])
 def test_mass_charge_energy(self):
  p=Parameters();r=predict(.6,parameters=p);mol=.6*p.effective_area_cells_cm2/(2*F)
  self.assertAlmostEqual(r['hydrogen_kg_h'],mol*.00201588*3600)
  self.assertAlmostEqual(r['oxygen_nm3_h']*2,r['hydrogen_nm3_h'])
  self.assertAlmostEqual(r['water_kg_h']/r['hydrogen_kg_h'],.01801528/.00201588)
  self.assertAlmostEqual(r['heat_kw']+mol*285830/1000,r['power_kw'])
 def test_nernst(self):self.assertGreater(predict(.5,pressure_bar=30)['voltage_v'],predict(.5,pressure_bar=1)['voltage_v'])
 def test_invalid(self):
  for j in [-1,float('nan'),float('inf')]:
   with self.assertRaises(ValueError):predict(j)
  with self.assertRaises(ValueError):predict(1,limiting_current=1)
 def test_thermal(self):
  self.assertAlmostEqual(thermal_step(30,10,0,ambient_c=30),30)
  self.assertAlmostEqual(thermal_step(20,10,1,ua_w_k=0,capacity_j_k=1000),30)
 def test_python_browser_parity(self):
  par=Parameters(**json.loads(Path('docs/validation.json').read_text())['parameters'])
  cases=[[0,20,1],[.6,55,30],[.83,50,34]]
  script="const m=require('./docs/model.js'); console.log(JSON.stringify("+json.dumps(cases)+".map(a=>m.predict(...a,"+json.dumps(vars(par))+"))));"
  result=json.loads(subprocess.check_output(['node','-e',script],text=True))
  for args,js in zip(cases,result):
   for k,v in predict(*args,parameters=par).items():
    if v is None:self.assertIsNone(js[k])
    else:self.assertAlmostEqual(v,js[k],places=10)
if __name__=='__main__':unittest.main()

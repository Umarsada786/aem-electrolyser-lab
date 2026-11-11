"""Reduced-order AEM model. SI except j [A/cm2], area [cm2], p [bar]."""
import math
from dataclasses import dataclass

F=96485.33212
R=8.314462618
VM=0.02241396954 # m3/mol, 0 C and 1 atm
@dataclass
class Parameters:
    activation_v: float=0.06
    exchange_a_cm2: float=0.001
    resistance_ohm_cm2: float=0.3
    temperature_coefficient_v_k: float=-0.002
    effective_area_cells_cm2: float=1463.0
    faradaic_efficiency: float=1.0

def predict(j,temperature_c=55,pressure_bar=30,parameters=None,oxygen_bar=1,auxiliary_kw=0,limiting_current=None):
    p=parameters or Parameters()
    values=[j,temperature_c,pressure_bar,oxygen_bar,auxiliary_kw,*vars(p).values()]
    if not all(math.isfinite(v) for v in values):raise ValueError('Inputs must be finite')
    if j<0 or not 0<=temperature_c<=100 or pressure_bar<=0 or oxygen_bar<=0 or auxiliary_kw<0:raise ValueError('Invalid operating conditions')
    if p.activation_v<0 or p.exchange_a_cm2<=0 or p.resistance_ohm_cm2<0 or p.effective_area_cells_cm2<=0 or not 0<p.faradaic_efficiency<=1:raise ValueError('Invalid parameters')
    if limiting_current is not None and (not math.isfinite(limiting_current) or limiting_current<=j):raise ValueError('Current must be below limiting current')
    t=temperature_c+273.15
    reversible=1.229-0.0009*(t-298.15)+R*t/(2*F)*math.log(pressure_bar*math.sqrt(oxygen_bar))
    activation=p.activation_v*(t/328.15)*math.asinh(j/(2*p.exchange_a_cm2))
    ohmic=p.resistance_ohm_cm2*j
    correction=p.temperature_coefficient_v_k*(temperature_c-55)*j
    transport=0 if limiting_current is None else -R*t/(2*F)*math.log(1-j/limiting_current)
    voltage=reversible+activation+ohmic+correction+transport
    if voltage<=0:raise ValueError('Parameters produce nonphysical voltage')
    ampere_cells=j*p.effective_area_cells_cm2
    power=ampere_cells*voltage/1000
    mol_s=ampere_cells*p.faradaic_efficiency/(2*F)
    kg_h=mol_s*0.00201588*3600
    return dict(voltage_v=voltage,power_kw=power,system_power_kw=power+auxiliary_kw,
      hydrogen_nm3_h=mol_s*VM*3600,hydrogen_kg_h=kg_h,oxygen_nm3_h=mol_s*VM*1800,
      water_kg_h=mol_s*0.01801528*3600,heat_kw=power-mol_s*285830/1000,
      efficiency_lhv=kg_h*33.33/(power+auxiliary_kw) if power+auxiliary_kw else None,
      specific_energy_kwh_kg=(power+auxiliary_kw)/kg_h if kg_h else None,
      reversible_v=reversible,activation_v=activation,ohmic_v=ohmic,temperature_correction_v=correction,transport_v=transport)

def thermal_step(temperature_c,dt_s,heat_kw,ambient_c=20,ua_w_k=10,capacity_j_k=50000,cooling_kw=0,heating_kw=0):
    """Exact constant-heat lumped balance; illustrative, NOT calibrated."""
    if not all(math.isfinite(x) for x in [temperature_c,dt_s,heat_kw,ambient_c,ua_w_k,capacity_j_k,cooling_kw,heating_kw]) or dt_s<0 or ua_w_k<0 or capacity_j_k<=0 or cooling_kw<0 or heating_kw<0:raise ValueError('Invalid thermal input')
    q=(heat_kw+heating_kw-cooling_kw)*1000
    if ua_w_k==0:return temperature_c+dt_s*q/capacity_j_k
    equilibrium=ambient_c+q/ua_w_k
    return equilibrium+(temperature_c-equilibrium)*math.exp(-ua_w_k*dt_s/capacity_j_k)

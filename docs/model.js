(function(root){
function predict(j,t,p,par,aux=0,limit=null){
 if(![j,t,p,aux,...Object.values(par)].every(Number.isFinite)||j<0||t<0||t>100||p<=0||aux<0||par.activation_v<0||par.exchange_a_cm2<=0||par.resistance_ohm_cm2<0||par.effective_area_cells_cm2<=0||par.faradaic_efficiency<=0||par.faradaic_efficiency>1)throw Error('Check operating conditions and parameters.');
 if(limit!==null&&(!Number.isFinite(limit)||limit<=j))throw Error('Current must be below limiting current.');
 const tk=t+273.15,F=96485.33212,R=8.314462618;
 const reversible=1.229-.0009*(tk-298.15)+R*tk/(2*F)*Math.log(p),activation=par.activation_v*tk/328.15*Math.asinh(j/(2*par.exchange_a_cm2)),ohmic=par.resistance_ohm_cm2*j,correction=par.temperature_coefficient_v_k*(t-55)*j,transport=limit===null?0:-R*tk/(2*F)*Math.log(1-j/limit),voltage=reversible+activation+ohmic+correction+transport;
 if(voltage<=0)throw Error('Parameters produce nonphysical voltage.');
 const ic=j*par.effective_area_cells_cm2,power=ic*voltage/1000,mol=ic*par.faradaic_efficiency/(2*F),kg=mol*.00201588*3600;
 return {voltage_v:voltage,power_kw:power,system_power_kw:power+aux,hydrogen_nm3_h:mol*.02241396954*3600,hydrogen_kg_h:kg,oxygen_nm3_h:mol*.02241396954*1800,water_kg_h:mol*.01801528*3600,heat_kw:power-mol*285830/1000,efficiency_lhv:power+aux?kg*33.33/(power+aux):null,specific_energy_kwh_kg:kg?(power+aux)/kg:null,reversible_v:reversible,activation_v:activation,ohmic_v:ohmic,temperature_correction_v:correction,transport_v:transport};
}
root.AEM={predict};if(typeof module!=='undefined')module.exports=root.AEM;
})(typeof window!=='undefined'?window:globalThis);

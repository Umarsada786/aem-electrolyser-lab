# AEM Lab — electrolyser modelling sandbox

A runnable, reduced-order anion exchange membrane (AEM) electrolyser model and a dependency-free browser simulator. Visitors can vary operating conditions, inspect polarization and energy outputs, download simulations and evaluate their own CSV measurements locally.

**Status:** calibrated and evaluated on a chronological holdout from one supplied experiment. This is not a fully validated multiphysics model, commercial digital twin or hardware controller.

## Try it

Open `docs/index.html` directly in any modern browser, or:

```bash
python -m http.server 8000 --directory docs
```

Visit http://localhost:8000. All calculations and CSV evaluation happen locally in the browser; no API keys or paid services are needed.

## Reproduce the validation

Requires Python 3.10+ and Node.js 18+ (Node is used for cross-language tests).

```bash
python -m venv .venv
source .venv/bin/activate # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python validate.py
python -m unittest discover -s tests -v
```

The script regenerates `docs/validation.json`, `data/predictions.csv` and `docs/dataset.js`. It always retains the train-only parameters used in the holdout evaluation. It never refits on the holdout. CI reruns the calibration and tests.

## Results

| Output | Calibration RMSE | Holdout RMSE | Holdout R² |
|---|---:|---:|---:|
| Cell voltage | 0.01594 V | 0.02979 V | 0.3426 |
| Stack power | 0.01032 kW | 0.02474 kW | 0.9940 |
| H₂ flow | 0.001117 Nm³/h | 0.001033 Nm³/h | 0.9998 |

386 source rows → one zero-current startup row excluded from electrochemical fitting → 269 calibration and 116 holdout observations. Calibration ends at 8,070 s; holdout begins at 8,100 s. No random split is used because adjacent observations are temporally correlated. All active startup observations are retained, including poorly fitted points.

These are conditional predictions given **measured current density, temperature and pressure**. They do not validate a free-running dynamic simulation. Holdout voltage R² is modest, despite small absolute error. Power and flow are strongly driven by current, and flow is nearly an algebraic multiple of current in the workbook; those metrics should not be interpreted as three independent confirmations.

The exchange-current fit reaches its lower bound (1e-8 A/cm²), and the fit Jacobian condition number is approximately 115,031. Coefficients are poorly identifiable and must not be interpreted as measured catalyst properties. There are no statistically calibrated prediction intervals.

## Model equations

Inputs: `j` in A/cm², `T` in K, hydrogen and oxygen pressures in bar, water activity assumed 1. The dataset pressure is assumed absolute; its reference is not specified in the workbook.

- `Erev = 1.229 - 0.0009*(T-298.15) + R*T/(2*F)*ln(pH2*sqrt(pO2))`
- `eta_act = a*(T/328.15)*asinh(j/(2*j0))`
- `eta_ohm = ASR*j`
- `eta_T = cT*(T-328.15)*j` (empirical correction, not a separate physical loss)
- `eta_transport = -R*T/(2*F)*ln(1-j/jlim)` (optional; disabled in calibration because no limiting current is identified)
- `Vcell = Erev + eta_act + eta_ohm + eta_T + eta_transport`
- `Pstack = j*(N*A)*Vcell/1000` kW
- `nH2 = etaF*j*(N*A)/(2*F)` mol/s; `nO2 = nH2/2`
- Net water consumption equals `nH2` mol/s, not circulation or make-up flow including losses.
- `Qnet = Pstack - nH2*285830/1000` kW (approximate 25 °C HHV energy balance; neglects sensible and latent flows).
- LHV efficiency uses 33.33 kWh/kg H₂ and includes the user-specified auxiliary electricity.

`N*A` is an effective aggregate cell-area product inferred using calibration power, voltage and current. Individual cell count and area cannot be recovered. Fitted product: 1463.20 cm². Faradaic efficiency inferred with 0 °C, 1 atm normal volume is 1.002585, which is nonphysical; it is capped at 1 and the resulting bias is retained. The mismatch may reflect rounding or a different normal-volume definition; it is not evidence of >100% efficiency.

At zero current the model returns equilibrium cell voltage and zero production/power. The workbook reports zero terminal voltage at shutdown, which is a different observable. Zero-current specific energy and efficiency are undefined where appropriate.

The Python `thermal_step` function integrates a constant-heat lumped thermal balance exactly over one time step. Heat capacity, UA, cooling and heating are user assumptions; thermal dynamics are **not fitted or validated**. No pressure-state dynamics, gas crossover, hydration, degradation, spatial fields or separate anode/cathode kinetics are resolved.

## Data provenance

Source: user-supplied **Supplementary data (2)(3).xlsx**, sheet `Exp`, rows 1–386. The first row includes both labels and a measurement; it is preserved. Column mapping and source SHA-256 are in `data/provenance.json`. Time is seconds, pressure converted MPa → bar by multiplication by 10. `data/experimental.csv` is the extracted experimental series with original worksheet row numbers.

The `Sim` and calculated comparison sheets are not experimental validation targets. `Energy balances` appears to represent a different aggregate scale, so it is not merged into the single-stack fit. Original workbook bytes are not bundled; extracted data are retained for reproducibility. The original publication/DOI and data reuse license were not supplied and are not invented. The code license does not grant rights to third-party experimental data.

Individual calibration ranges: j = 0.06143–0.83583 A/cm², temperature = 20.703–55.115 °C, pressure = 1–33.964 bar. These ranges are a bounding box, not proof that every combination was tested. Temperature, pressure and current are correlated in one trajectory. Extrapolation and altered parameters are marked in the UI.

## Python usage

```python
import json
from model import Parameters, predict
p = Parameters(**json.load(open('docs/validation.json'))['parameters'])
result = predict(0.6, temperature_c=55, pressure_bar=30, parameters=p)
print(result)
```

## Publish on GitHub Pages

After creating a public GitHub repository and pushing these files, go to **Settings → Pages → Build and deployment → Deploy from a branch**, choose **main** and **/docs**, and Save. GitHub will display the live URL after deployment. No build workflow or secrets are necessary for the website. The test workflow is separate from Pages deployment.

For an authenticated GitHub CLI:

```bash
git init -b main
git add .
git commit -m "Add AEM model, experimental validation and browser simulator"
gh repo create aem-electrolyser-lab --public --source=. --remote=origin --push
gh api --method POST repos/$(gh api user --jq .login)/aem-electrolyser-lab/pages -f 'source[branch]=main' -f 'source[path]=/docs'
```

If Pages is already configured, use its Settings page rather than creating it again.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Useful next evidence: independently measured steady-state polarization curves across temperature/pressure, EIS resistance, active area and cell count, gas reference conditions, Faradaic efficiency and replicate experiments. These are necessary to move beyond effective coefficients.

## Scientific background

Chang, S.-C., Gu, R.-E. & Chan, Y.-H. (2024). *Parameter Analysis of Anion Exchange Membrane Water Electrolysis System by Numerical Simulation*. Energies 17(22), 5682. https://doi.org/10.3390/en17225682

The paper motivates the polarization decomposition and system perspective. This repository implements a simpler effective model with its own calibration; it does not reproduce the paper's full model or claim its validation results.

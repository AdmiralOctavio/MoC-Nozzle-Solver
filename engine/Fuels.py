from dataclasses import dataclass, field
from typing import Optional
import numpy as np
from rocketcea.cea_obj_w_units import CEA_Obj
from rocketcea.blends import newFuelBlend
import cea


"""
Fuel definitions
Feel free to create your own fuel blends depending on what you need to simulate
Parameters.py utilises the new and old cea libraries, so both inputs must be given

"""


@dataclass
class PropellantConfig:
    # For new CEA
    name: str
    reac_names: list
    fuel_weights: np.ndarray
    oxidant_weights: np.ndarray
    T_reactant: np.ndarray

    # For old rocketcea
    ox_name: str = ""
    fuel_name: str = ""

    def make_mixture(self):
        # Return (reac, prod, weights_fn) for new cea
        reac = cea.Mixture(self.reac_names)
        prod = cea.Mixture(self.reac_names, products_from_reactants=True)
        return reac, prod

    def get_weights(self, of_ratio: float) -> np.ndarray:
        # Convert O/F ratio to mass fractions using this propellant's weight vectors.
        reac, _ = self.make_mixture()
        return reac.of_ratio_to_weights(
            self.oxidant_weights, self.fuel_weights, of_ratio
        )

    def get_enthalpy_coeff(self, of_ratio: float) -> float:
        """Return hc (enthalpy / R) for this propellant at the given O/F."""
        reac, _ = self.make_mixture()
        weights = self.get_weights(of_ratio)
        return reac.calc_property(cea.ENTHALPY, weights, self.T_reactant) / cea.R

    def make_cea_obj(self) -> CEA_Obj:
        """Return a rocketcea CEA_Obj configured for this propellant."""
        return CEA_Obj(
            oxName=self.ox_name,
            fuelName=self.fuel_name,
            temperature_units="degK",
            pressure_units="bar",
        )


_ethanol100 = newFuelBlend(["Ethanol", "H2O"], [100, 0])
_ethanol80 = newFuelBlend(["Ethanol", "H2O"], [80, 20])

N2O_Ethanol80 = PropellantConfig(
    name="N2O / Ethanol (80%)",
    reac_names=[b"C2H5OH(L)", b"H2O", b"N2O"],
    fuel_weights=np.array([0.8, 0.2, 0.0]),
    oxidant_weights=np.array([0.0, 0.0, 1.0]),
    T_reactant=np.array([298.15, 298.15, 298.15]),
    ox_name="N2O",
    fuel_name=_ethanol80,
)

N2O_Ethanol100 = PropellantConfig(
    name="N2O / Ethanol (100%)",
    reac_names=[b"C2H5OH(L)", b"H2O", b"N2O"],
    fuel_weights=np.array([1, 0, 0.0]),
    oxidant_weights=np.array([0.0, 0.0, 1.0]),
    T_reactant=np.array([298.15, 298.15, 298.15]),
    ox_name="N2O",
    fuel_name=_ethanol100,
)

LOX_RP1 = PropellantConfig(
    name="LOX / RP-1",
    reac_names=[b"RP-1", b"O2(L)"],
    fuel_weights=np.array([1.0, 0.0]),
    oxidant_weights=np.array([0.0, 1.0]),
    T_reactant=np.array([298.15, 90.15]),
    ox_name="LOX",
    fuel_name="RP-1",
)

LOX_LH2 = PropellantConfig(
    name="LOX / LH2",
    reac_names=[b"H2(L)", b"O2(L)"],
    fuel_weights=np.array([1.0, 0.0]),
    oxidant_weights=np.array([0.0, 1.0]),
    T_reactant=np.array([20.15, 90.15]),
    ox_name="LOX",
    fuel_name="LH2",
)

LOX_Ethanol100 = PropellantConfig(
    name        = "LOX / Ethanol (100%)",
    reac_names  = [b"C2H5OH(L)", b"O2(L)"],
    fuel_weights    = np.array([1.0, 0.0]),   
    oxidant_weights = np.array([0.0, 1.0]),   
    T_reactant  = np.array([298.15, 90.15]),  
    ox_name     = "LOX",
    fuel_name   = _ethanol100,
)

LOX_Ethanol80 = PropellantConfig(
    name        = "LOX / Ethanol (80%)",
    reac_names  = [b"C2H5OH(L)", b"H2O(L)", b"O2(L)"],  
    fuel_weights    = np.array([0.8, 0.2, 0.0]),   
    oxidant_weights = np.array([0.0, 0.0, 1.0]),    
    T_reactant  = np.array([298.15, 298.15, 90.15]),
    ox_name     = "LOX",
    fuel_name   = _ethanol80,
)


PROPELLANTS: dict[str, PropellantConfig] = {
    "n2o_ethanol100": N2O_Ethanol100,
    "n2o_ethanol80": N2O_Ethanol80,
    "lox_rp1": LOX_RP1,
    "lox_lh2": LOX_LH2,
    "lox_ethanol100": LOX_Ethanol100,
    "lox_ethanol80": LOX_Ethanol80,
}

PROPELLANT_OPTIONS = [(k, v.name) for k, v in PROPELLANTS.items()]

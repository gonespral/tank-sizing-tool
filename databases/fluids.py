import valispace

fluids = {
    "lox": {
        "name": "Liquid Oxygen",
        "density": valispace.get_vali(13371)["value"],
        "thermal_conductivity": 0.026,
        "specific_heat": 2000,
        "thermal_expansion_coeff": 0.00004
    },
    "ethanol": {
        "name": "Ethanol",
        "density": valispace.get_vali(13365)["value"],
        "thermal_conductivity": 0.24,
        "specific_heat": 2500,
        "thermal_expansion_coeff": 0.00004
    }
}

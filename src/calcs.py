from mace.calculators import MACECalculator


def calc():

    mlcalculator = MACECalculator(
        model_paths='/global/home/users/kumaranu/Documents/gpu_jobs/MACE_model.model',
        device='cuda',
        default_dtype="float64",
    )
    return mlcalculator

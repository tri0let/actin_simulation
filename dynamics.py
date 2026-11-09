import random
import math
import matplotlib.pyplot as plt

def sample_exp(mean: float) -> float:
    return -math.log(1 - random.random()) * mean


from enum import Enum
from dynamics import *
from typing import Literal

class ActinState(Enum):
    adp = '.'
    atp = 'O'

class Actin:
    lifetime = 0.5

    def __init__(self, state: int) -> None:
        if state:
            self.state = ActinState.atp
            self.time = sample_exp(Actin.lifetime)
        else:
            self.state = ActinState.adp
            self.time = -1.0

    def hydrolyze(self):
        self.state = ActinState.adp

    def update(self, dt: float):
        if self.state == ActinState.atp:
            self.time -= dt
            if self.time <= 0:
                self.hydrolyze()

    def __str__(self):
        return self.state.value

class ActinPolymer:

    Cc = {
        ActinState.atp: 1.0,
        ActinState.adp: 2.0
    }    
    
    K: dict[str, dict[ActinState, dict[str, float]]] = {
        'plus': {
            ActinState.atp: {'on': 12.0},
            ActinState.adp: {'on': 4.0},
        },
        'minus': {
            ActinState.atp: {'on': 1.0},
            ActinState.adp: {'on': 2.0},
        },
    }

    for val in K.values():
        for state in val.keys():
            val[state]['off'] = val[state]['on'] * Cc[state]

    def __init__(self, monomers: list[Actin], C_monomer: float) -> None:
        self.monomers = monomers
        self.C_monomer = C_monomer
        self.times = {
            key: {
                state: {
                    k: sample_exp(1 / (v * self.C_monomer)) if k == 'on' else sample_exp(1 / v) for k, v in ActinPolymer.K[key][state].items()
                } for state in ActinPolymer.K[key].keys()
            } for key in ActinPolymer.K.keys()
        }

        self.added = {'plus': 0, 'minus': 0}

    def reset_times(self, end: str, state: ActinState, dir: str):
        if end not in ['plus', 'minus'] or dir not in ['on', 'off']:
            raise ValueError
        if dir == 'on':
            self.times[end][state][dir] = sample_exp(1 / (ActinPolymer.K[end][state][dir] * self.C_monomer))
        else:
            self.times[end][state][dir] = sample_exp(1 / (ActinPolymer.K[end][state][dir]))

    def add(self, monomer: Actin, end: str):
        if end == 'minus':
            self.monomers.insert(0, monomer)
        elif end == 'plus':
            self.monomers.append(monomer)
        else: raise ValueError
        self.added[end] += 1


    def remove(self, end: str):
        if len(self.monomers) == 1:
            return
        if end == 'minus':
            self.monomers.pop(0)
        elif end == 'plus':
            self.monomers.pop(-1)
        else: raise ValueError
        self.added[end] -= 1

    def end_state(self, end: str) -> ActinState:
        if end == 'plus':
            return self.monomers[-1].state
        elif end == 'minus':
            return self.monomers[0].state
        else: raise ValueError

    def second_state(self, end:str) -> ActinState:
        if len(self.monomers) < 2:
            return self.end_state(end)
        if end == 'plus':
            return self.monomers[-2].state
        elif end == 'minus':
            return self.monomers[1].state
        else: raise ValueError

    def update(self, dt: float):
        for actin in self.monomers:
            actin.update(dt)

        for end in ['plus', 'minus']:
            for state in ActinState:
                for dir in ['on', 'off']:
                    self.times[end][state][dir] -= dt
                    if self.times[end][state][dir] <= 0:
                        self.reset_times(end, state, dir)
                        if self.end_state(end) == state and dir == 'on':
                            self.add(Actin(1), end)
                        if self.second_state(end) == state and dir == 'off':
                            if end == 'minus' and state == ActinState.adp:
                                pass
                            self.remove(end)


    def __str__(self):
        return '-' + ''.join([str(m) for m in self.monomers]) + '+'

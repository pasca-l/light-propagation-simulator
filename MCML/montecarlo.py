import math
import numpy as np
import random


class Calculator:
    def randnum(self):
        return random.random()

    def sign(self, num):
        return (num > 0) - (num < 0)

    def step_size(self, mu_a, mu_s):
        return -1 * math.log(self.randnum()) / (mu_a + mu_s)

    def delta_weight(self, mu_a, mu_s, weight):
        return mu_s * weight / (mu_a + mu_s)

    # calculation differ from paper
    def border_dist(self, z, uz, z_upper, z_lower):
        if uz < 0:
            return (z - z_upper)
        elif uz > 0:
            return (z_lower - z)
        else:
            return float('inf')

    def ref_specular(self, n1, n0=1):
        return ((n0 - n1)*(n0 - n1)) / ((n0 + n1)*(n0 + n1))

    def ref_at_border(self, direction, n1, n2):
        ux, uy, uz = direction.values()

        if abs(uz) == 1:
            a_i = 2 * math.pi
            a_t = 2 * math.pi
        else:
            a_i = math.acos(abs(uz))
            a_t = math.asin(n1 / n2 * math.sin(a_i))

        ref = 0.5 * ((math.sin(a_i - a_t) * math.sin(a_i - a_t)
                    / math.sin(a_i + a_t) * math.sin(a_i + a_t))
                    + (math.tan(a_i - a_t) * math.tan(a_i - a_t)
                    / math.tan(a_i + a_t) * math.tan(a_i + a_t)))

        if self.randnum() > ref:
            ux = math.sin(a_t) * ux / math.sin(a_i)
            uy = math.sin(a_t) * uy / math.sin(a_i)
            uz = self.sign(uz) * math.cos(a_t)
            return {'x':ux, 'y':uy, 'z':uz}

        else:
            return {'x':ux, 'y':uy, 'z':-1*uz}

    def scat_dir(self, direction, g=0):
        ux, uy, uz = direction.values()

        if g == 0:
            theta = math.acos(2 * self.randnum() - 1)
        phi = 2 * math.pi * self.randnum()

        if abs(uz) > 0.99999:
            ux = math.sin(theta) * math.cos(phi)
            uy = math.sin(theta) * math.cos(phi)
            uz = self.sign(uz) * math.cos(theta)
            return {'x':ux, 'y':uy, 'z':uz}

        else:
            ux = (math.sin(theta)
                * (ux * uz * math.cos(phi) - uy * math.sin(phi))
                / math.sqrt(1 - uz * uz)
                + ux * math.cos(theta))
            uy = (math.sin(theta)
                * (uy * uz * math.cos(phi) - ux * math.sin(phi))
                / math.sqrt(1 - uz * uz)
                + uy * math.cos(theta))
            uz = (-1 * math.sin(theta) * math.cos(phi) * math.sqrt(1 - uz * uz)
                + uz * math.cos(theta))
            return {'x':ux, 'y':uy, 'z':uz}

    def russian_roulette(self, weight, chance=10):
        if self.randnum() <= 1 / chance:
            return chance * weight
        else:
            return 0


class Model:
    def __init__(self):
        self.model = self.load_model()

    def load_model(self):
        # data in [[mu_a (mm-1), mu_s (mm-1), n, z_upper, z_lower], ...]
        return [[0.02, 10.0, 1.4, 0, 1],
                [0.018, 1.9, 1.1, 1, 2],
                [0.018, 1.9, 1.3, 2, 3]]

    def optic_param(self, z, uz, mu_a=0, mu_s=1, n=2, z_up=3, z_low=4):
        if uz > 0:
            for i, layer in enumerate(self.model):
            # if z coordinate within upper and lower depth
                if layer[z_up] <= z and layer[z_low] > z:
                    return {'mu_a':layer[mu_a],
                            'mu_s':layer[mu_s],
                            'n':layer[n],
                            'z_up':layer[z_up],
                            'z_low':layer[z_low],
                            'index':i}
        elif uz < 0:
            for i, layer in enumerate(self.model):
                if layer[z_up] < z and layer[z_low] >= z:
                    return {'mu_a':layer[mu_a],
                            'mu_s':layer[mu_s],
                            'n':layer[n],
                            'z_up':layer[z_up],
                            'z_low':layer[z_low],
                            'index':i}

    def model_border(self, z_up=3, z_low=4):
        z_min, z_max = 0, 0
        for layer in self.model:
            if layer[z_up] < z_min:
                z_min = layer[z_up]
            if layer[z_low] > z_max:
                z_max = layer[z_low]
        return z_min, z_max


class Recorder:
    def __init__(self):
        self.absorbance = np.empty(0)
        self.reflectance = np.empty(0)
        self.transmittance = np.empty(0)
        self.path = np.empty(0)

    def record(self, *data, abs=False, ref=False, tra=False):
        if abs:
            self.absorbance = np.append(self.absorbance, data)
        if ref:
            self.reflectance = np.append(self.reflectance, data)
        if tra:
            self.transmittance = np.append(self.transmittance, data)

    def optimize(self):
        # method to sum up values in close range
        return


class Simulator:
    def __init__(self, model, calculator, recorder):
        self.model = model
        self.calc = calculator
        self.recorder = recorder
        self.coord = {'x':0, 'y':0, 'z':0}
        self.direc = {'x':0, 'y':0, 'z':1}
        self.weight = 1
        self.w_thresh = 0.0001
        self.optic_param = None
        self.dead = False

    def show_step(self):
        print("photon at:    ", self.coord)
        print("directing to: ", self.direc)
        print("with weight:  ", self.weight)
        print("in medium number: " + str(self.optic_param['index']))
        print("")

    def light_propagation(self):
        step = 0

        self.optic_param = self.model.optic_param(self.coord['z'],
                                                  self.direc['z'])
        self.weight = 1 - self.calc.ref_specular(self.optic_param['n'])

        while not self.dead:
            if step == 0:
                step = self.calc.step_size(self.optic_param['mu_a'],
                                           self.optic_param['mu_s'])
            border = self.calc.border_dist(self.coord['z'],
                                           self.direc['z'],
                                           self.optic_param['z_up'],
                                           self.optic_param['z_low'])

            # if does not hit boundary
            if border > step:
                # update parameters
                for k in self.coord.keys():
                    self.coord[k] = self.coord[k] + self.direc[k] * step
                self.weight = (self.weight - self.calc.delta_weight(
                                                self.optic_param['mu_a'],
                                                self.optic_param['mu_s'],
                                                self.weight))
                self.direc.update(self.calc.scat_dir(self.direc))
                step = 0

                # record absorbance
                self.recorder.record(
                    self.coord,
                    self.calc.delta_weight(self.optic_param['mu_a'],
                                           self.optic_param['mu_s'],
                                           self.weight),
                    abs=True)


            elif border <= step:
                # update parameters
                for k in self.coord.keys():
                    self.coord[k] = (self.coord[k]
                                    + border * self.calc.sign(self.direc[k]))
                # ref_at_border must take in next layer refractive index
                self.direc.update(
                    self.calc.ref_at_border(self.direc,
                                            self.optic_param['n'],
                                            1.4))
                # below step must change to consider different optical parameters
                step = step - border
                # -------------
                self.optic_param = self.model.optic_param(self.coord['z'],
                                                          self.direc['z'])

                # if photon directs outside of model
                if self.optic_param == None:
                    z_min, z_max = self.model.model_border()
                    if self.coord['z'] == z_min:
                        # record reflectance
                        self.recorder.record(self.coord, self.weight, ref=True)
                    if self.coord['z'] == z_max:
                        # record transmittance
                        self.recorder.record(self.coord, self.weight, tra=True)

                    self.dead = True
                    break

            # weight check
            if self.weight < self.w_thresh:
                self.weight = self.calc.russian_roulette(self.weight)

            if self.weight == 0:
                self.dead = True
                break


if __name__ == '__main__':
    model = Model()
    calculator = Calculator()
    recorder = Recorder()

    count = 5
    while count > 0:
        simulator = Simulator(model, calculator, recorder)
        simulator.light_propagation()
        count -= 1

    print(simulator.recorder.reflectance)

#!/usr/bin/python
# -*- coding: ISO-8859-1 -*-
"""
Mara.ANN.Learn (v. 0.1):
    Learning algorithms for neural networks.

    Author: Martti Louhivuori (martti.louhivuori@helsinki.fi)

    Date: 14.2.2006
"""

class Teacher:

    def __init__(self, net, method='bp', targets=None, rate=1,
                gradient_delta=0.1):
        self.net = net
        self.method = method
        self.targets = targets
        self.rate = rate
        self.__valid_methods__ = ['bp']
        self.__gradient_delta__ = gradient_delta

    def get_valid_methods(self):
        return self.__valid_methods__

    def set_gradient_delta(self, delta):
        self.__gradient_delta__ = delta

    def get_gradient_delta(self):
        return self.__gradient_delta__

    def teach(self, net=None, targets=None):
        if net is None:
            net = self.net
        if targets is None:
            targets = self.targets
        errors = []
        for neuron,target in zip(net.layers['output'], targets):
            errors.append(float(target-neuron.status()))
        total_error = 0
        for e in errors:
            total_error += e**2
        total_error /= 2
        done = []
        current = []
        for neuron,error in zip(net.layers['output'],errors):
            signal = neuron.input_signal()
            neuron.gradient = error*neuron.evolve(signal+self.__gradient_delta__)/self.__gradient_delta__
            done.append(neuron)
            for k in neuron.inputs.keys():
                if k not in current:
                    current.append(k)
        while len(current) > 0:
            next = []
            for n in current:
                if n not in done:
                    signal = n.input_signal()
                    delta = n.evolve(signal+self.__gradient_delta__)/self.__gradient_delta__
                    sum = 0
                    for o in n.outputs:
                        sum += o.gradient*o.inputs[n]['weight']
                        o.inputs[n]['weight'] += self.rate*o.gradient*n.status()
                    n.gradient = delta*sum
                    for k in n.inputs.keys():
                        if k not in next:
                            next.append(k)
                    done.append(n)
            current = next

            
            tree = []
            for 



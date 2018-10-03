#!/usr/bin/python
# -*- coding: ISO-8859-1 -*-
"""
Mara.ANN.Arch (v. 0.1):
    Neural network architectures.

    Author: Martti Louhivuori (martti.louhivuori@helsinki.fi)

    Date: 14.2.2006
"""
from random import Random

class Neuron:
    "A single artificial neuron."

    def __init__(self, id=None, activation='linear', 
                negative_weights=False, max_history_size=10, 
                random=Random(), bias=False):
        """
        Creates a neuron of given type.
        """
        self.__negative_weights__ = negative_weights
        self.__random__ = random
        self.__max_history_size__ = max_history_size
        self.id = id
        self.history = [1]
        self.inputs = {}
        self.outputs = []
        self.set_activation_fn(activation)
        self.bias = False
        self.gradient = None
        if bias:
            self.__max_history_size__ = 1
            self.bias = True
    
    def set_activation_fn(self, fn):
        """
        Sets a new activation function for the neuron.
        """
        functions = ['linear', 'step', 'sigmoid']
        if functions.count(fn):
            self.__activation__ = fn
            return True
        else:
            print "Unknown activation funtion: " + str(fn)
            return False

    def get_activation_fn(self):
        return self.__activation__

    def set_negative_weights(self, value):
        if value:
            self.__negative_weights__ = 1
        else:
            self.__negative_weights__ = 0

    def get_negative_weights(self):
        return self.__negative_weights__

    def add_input(self, other, weight=None, delay=0, fixed=False):
        """
        Adds an input connection from another neuron.
        """
        if self.bias:
            print "Can't add a connection to a bias node."
            return False
        else:
            if weight is None:
                if self.__negative_weights__:
                    weight = self.__random__.uniform(-1,1)
                else:
                    weight = self.__random__.uniform(0,1)
            if delay < 0:
                print "Delay can be only positive or zero."
                return False
            if self.inputs.has_key(other):
                print "There is already a connection between " + str(self.id) + " and " + str(other.id) + ". Only adjusting parameters."
            self.inputs[other] = {'weight':weight, 'delay':delay, 'fixed':fixed}
            return True

    def add_output(self, other):
        """
        Adds an output connection to another neuron.
        """
        if self.outputs.count(other):
            print "There is already a connection to " + str(other.id)
        else:
            self.outputs.append(other)

    def is_connected(self, other):
        """
        Tells what kind of a connection, if any, there is to an other neuron.
        """
        if self.inputs.has_key(other):
            return 1
        elif self.outputs.count(other):
            return -1
        else:
            return 0

    def status(self, delay=0):
        if delay < len(self.history):
            return self.history[-1*(delay+1)]
        else:
            print "Not long enough history for neuron " + str(self.id) + "."
            return 0

    def input_signal(self):
        sum = 0
        for neuron,link in self.inputs.items():
            sum += neuron.status(link['delay'])*link['weight']
        return sum

    def evolve(self, signal=None):
        """
        Calculate a new state for the neuron.
        """
        if signal is None:
            if len(self.history) > self.__max_history_size__ - 1:
                self.history.pop(0)
            if self.bias:
                self.history.append(1)
                return True
        if self.bias:
            return 1
        else:
            sum = self.input_signal()
            if sum != 0:
                if self.__activation__ is 'linear':
                    new = sum/len(self.inputs)
                elif self.__activation__ is 'step':
                    if sum > 0:
                        new = 1
                    else:
                        new = 0
                elif self.__activation__ is 'sigmoid':
# KORJAA!!!
                    new = sum/len(self.inputs)
                else:
                    print "Activation function " + str(self.__activation__) + " not defined."
                    new = None
            else:
                return True
            if signal is None:
                self.history.append(new)
                if new is None:
                    return False
                else:
                    return True
            else:
                return new

    def adjust_connection(self, other, weight=None, delay=None, fixed=None):
        """
        Adjusts the parameters of a connection.
        """
        if self.inputs.has_key(other):
            link = self.inputs[other]
            if weight is not None:
                link['weight'] = weight
            if delay is not None:
                if delay < 0:
                    print "Delay can only be a positive integer."
                else:
                    link['delay'] = delay
            if fixed is not None:
                if fixed:
                    link['fixed'] = True
                else:
                    link['fixed'] = False
            self.inputs[other] = link
            return True
        else:
            print "Unknown connection: " + str(other.id)

class Network:
    "Provides all the functionality needed to create and upkeep a neural network."

    def __init__(self, seed=None, defaults=None, wmin=0, wmax=1):
        if defaults is None:
            self.__defaults__ = {'activation':'linear',
                                'negative_weights':False,
                                'max_history_size':10}
        else:
            # Pitäisi varmaan tarkistaa rakenne...
            self.__defaults__ = defaults
        self.__random__ = Random(seed)
        self.__weight_min__ = wmin
        self.__weight_max__ = wmax
        self.neurons = []
        self.layers = {'input' : [], 'output' : []}
        self.biases = []
        
    def add_neuron(self, layer, neuron=None, bias=False):
        if neuron is None:
            neuron = Neuron(id = len(self.neurons)+1,
                            activation = self.__defaults__['activation'],
                            negative_weights = self.__defaults__['negative_weights'],
                            max_history_size = self.__defaults__['max_history_size'],
                            bias = bias)
        if self.layers.has_key(layer):
            self.layers[layer].append(neuron)
        else:
            self.layers[layer] = [neuron]
        self.neurons.append(neuron)
        if bias:
            self.biases.append(neuron)

    def remove_neuron(self, layer, neuron):
        self.layer[layer].remove(neuron)

    def connect_layers(self, first, second, two_way=False,
                        weights=None, delays=None):
        if weights is not None:
            weight = weights
        else:
            if self.__defaults__['negative_weights']:
                weight = self.__random__.uniform(-1,1)
            else:
                weight = self.__random__.uniform(0,1)
        for sender in self.layers[first]:
            for receiver in self.layers[second]:
                receiver.add_input(sender)
                sender.add_output(receiver)
                if two_way:
                    sender.add_input(receiver)
                    receiver.add_output(sender)

    def status(self):
        max = 0
        names = ['input']
        values = [[]]
        for k,v in self.layers.items():
            a = []
            for n in v:
                s = str(n.status())
                if n.bias:
                    s += "b"
                a.append(s)
            if k is 'input':
                values[0] = a
            elif k is 'output':
                o = a
            else:
                names.append(k)
                values.append(a)
            if len(v) > max:
                max = len(v)
        names.append('output')
        values.append(o)
        for i in range(len(values)):
            while len(values[i]) < max:
                values[i].append("")
        s = (max+1)*[""]
        for n,v in zip(names,values):
            s[0] += str(n)+"\t"
            for i,j in zip(v,range(1,len(v)+1)):
                s[j] += str(i)+"\t"
        for r in s:
            r.rstrip("\t")
            print r


class FeedForwardNetwork(Network):
    "A feed-forward network."

    def __init__(self, input = None, output = None, layers = None, 
                size = None):
        Network.__init__(self)
        i = 2
        o = 1
        # Check input and parse layer sizes.
        if input is not None and type(input) is type(1) and input > 0:
            i = input
        else:
            print "Invalid value given as the number of input layers."
        if output is not None and type(output) is type(1) and output > 0:
            o = output
        else:
            print "Invalid value given as the number of output layers."
        m = {}
        if layers is not None:
            if type(layers) is type({}):
                for n,v in layers.items():
                    if type(v) is type(1) and v > 0:
                        m[n] = v
                    else:
                        print "Invalid value given as node count in layer definitions."
            else:
                print "Invalid data format used for layers."
        if size is not None and type(size) is type([]) and len(size) > 1:
            if type(size[0]) is type(1) and size[0] > 0:
                i = size[0]
            else:
                print "Invalid value given as node count in size definitions."
            if type(size[-1]) is type(1) and size[-1] > 0:
                o = size[-1]
            else:
                print "Invalid value given as node count in size definitions."
            m = {}
            for n,v in zip(range(1,len(size)-1),size[1:-1]):
                if type(v) is type(1) and v > 0:
                    m[str(n)] = v
        # Create layers.
        names = ['input']
        values = [i]
        for n,v in m.items():
            names.append(n)
            values.append(v)
        names.append('output')
        values.append(o)
        for n,v in zip(names,values):
            while len(self.layers[n]) < v:
                self.add_neuron(n)
            while len(self.layers[n]) > v:
                self.remove_neuron(n)
        # Add biases.
        for n in self.layers.keys():
            if n is not 'output':
                self.add_neuron(n,bias=True)
        # Connect layers.
        for i in range(len(names)-1):
            self.connect_layers(names[i],names[i+1])


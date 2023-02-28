#!/usr/bin/env python
# coding: utf-8

# In[81]:


import cirq
from math import pi
import numpy as np


# In[82]:


#pair generation circuit class

class pair_generation(cirq.Gate):
    def __init__(self, theta):
        super (pair_generation, self)
        self.theta = theta
        
    def _num_qubits_(self):
        return 4
    
    def _decompose_(self, qubits):
        q0, q1, q2, q3 = qubits
        theta = self.theta
        
        yield cirq.CX(q1, q0)
        yield cirq.CX(q3, q1)
        yield cirq.CZ(q0, q1)
        yield cirq.H(q0)
        yield cirq.CX(q3, q2)
        yield cirq.ry(rads= theta).on(q2)
        yield cirq.ry(rads=- theta).on(q3)
        yield cirq.CX(q0, q3)
        yield cirq.CX(q0, q2)
        yield cirq.ry(rads=- theta).on(q3)
        yield cirq.ry(rads= theta).on(q2)
        yield cirq.CX(q1, q2)
        yield cirq.CX(q2, q3)
        yield cirq.ry(rads=- theta).on(q2)
        yield cirq.ry(rads= theta).on(q3)
        yield cirq.CX(q0, q3)
        yield cirq.CZ(q0, q3)
        yield cirq.ry(rads=- theta).on(q2)
        yield cirq.ry(rads= theta).on(q3)
        yield cirq.CX(q3, q2)
        yield cirq.CX(q1, q3)
        yield cirq.H(q3)
        yield cirq.CX(q3, q1)
        yield cirq.CX(q1, q0)
        
    def _circuit_diagram_info_(self, args):
        return ["Px"]*self._num_qubits_()

        


# In[83]:

# In[84]:


#given rotation generation

class SqrtISWAP(cirq.Gate):
    def __init__(self):
        super(SqrtISWAP, self)

    def _num_qubits_(self):
        return 2
    
    # Obtained from: https://github.com/quantumlib/Cirq/blob/a22269dfe41b0da78243bbd210a915d26cc7d25f/cirq-core/cirq/ops/swap_gates.py#L166
    def _decompose_(self, qubits):
        q0, q1 = qubits
        yield cirq.CX(q0,q1)
        yield cirq.H(q0)
        yield cirq.CX(q1,q0)
        yield cirq.ZPowGate(exponent=0.25).on(q0)
        yield cirq.CX(q1,q0)
        yield cirq.ZPowGate(exponent=-0.25).on(q0)
        yield cirq.H(q0)
        yield cirq.CX(q0,q1)
        
    def _circuit_diagram_info_(self, args):
        return ["sqrt(iSWAP)"] * self.num_qubits()


# In[85]:


class IonSupportedGivens(cirq.Gate):
    def __init__(self, theta):
        super(IonSupportedGivens, self)
        self.theta = theta

    def _num_qubits_(self):
        return 2

    def _decompose_(self, qubits):
        
        # Have to EXPLICITLY call decomposition, Azure Quantum doesn't seem to auto-decompose
        # to IonQ supported gates
        q0, q1 = qubits
        yield SqrtISWAP().on(q0,q1)._decompose_()
        
        yield cirq.rz(rads= -self.theta).on(q0)
        yield cirq.rz(rads=self.theta + np.pi).on(q1)

        yield SqrtISWAP().on(q0,q1)._decompose_()
        yield cirq.rz(rads=np.pi).on(q1)
        
    def _circuit_diagram_info_(self, args):
        return ["IonGivens({:.3f})".format(self.theta)] * self.num_qubits()


# In[86]:


#make a circuit with a random angle to test

#https://pubs.rsc.org/en/content/articlelanding/2022/sc/d1sc05691c#cit39
#from fig 1
class vqe(cirq.Gate):
    def __init__(self, theta):
        super(vqe, self)
        self.theta = theta
        
    def _num_qubits_(self):
        return 8
    
    def _decompose_(self, qubits):
        q0, q1, q2, q3, q4, q5, q6, q7 = qubits
        
        yield IonSupportedGivens(theta = self.theta).on(q0, q2)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q4, q6)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q1, q3)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q5, q7)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q2, q4)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q3, q5)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q0, q2)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q4, q6)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q1, q3)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q5, q7)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q2, q4)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q3, q5)._decompose_()
        yield pair_generation(theta = self.theta).on(q0, q1, q2, q3)._decompose_()
        yield pair_generation(theta = self.theta).on(q4, q5, q6, q7)._decompose_()
        yield pair_generation(theta = self.theta).on(q2, q3, q4, q5)._decompose_()
        
        yield pair_generation(theta = self.theta).on(q0, q1, q2, q3)._decompose_()
        yield pair_generation(theta = self.theta).on(q4, q5, q6, q7)._decompose_()
        yield pair_generation(theta = self.theta).on(q2, q3, q4, q5)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q0, q2)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q4, q6)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q1, q3)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q5, q7)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q2, q4)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q3, q5)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q0, q2)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q4, q6)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q1, q3)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q5, q7)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q2, q4)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q3, q5)._decompose_()
        
        
    def _circuit_diagram_info_(self, args):
        return ["VQE".format(self.theta)] * self.num_qubits()


class fourQubitVQE(cirq.Gate):
    def __init__(self, theta):
        super(fourQubitVQE, self)
        self.theta = theta
        
    def _num_qubits_(self):
        return 4
    
    def _decompose_(self, qubits):
        q0, q1, q2, q3 = qubits
        
        yield IonSupportedGivens(theta = self.theta).on(q0, q2)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q1, q3)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q0, q2)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q1, q3)._decompose_()
        yield pair_generation(theta = self.theta).on(q0, q1, q2, q3)._decompose_()
        
        yield pair_generation(theta = self.theta).on(q0, q1, q2, q3)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q0, q2)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q1, q3)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q0, q2)._decompose_()
        yield IonSupportedGivens(theta = self.theta).on(q1, q3)._decompose_()
        
        
    def _circuit_diagram_info_(self, args):
        return ["FourQubitVQE".format(self.theta)] * self.num_qubits()
        


# In[87]:


# test circuit 

if __name__ == "__main__":
    angle = np.pi/8

    vqe_cir = cirq.Circuit()
    q0, q1, q2, q3 = cirq.LineQubit.range(4)

    vqe_cir.append(fourQubitVQE(theta = angle).on(q0, q1, q2, q3)._decompose_())
    print(vqe_cir)


# In[ ]:





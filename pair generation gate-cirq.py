#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cirq
from math import pi


# In[2]:


circuit= cirq.Circuit()
q0, q1, q2, q3= cirq.LineQubit.range(4)

circuit.append(cirq.CNOT(q1, q0))
circuit.append(cirq.CNOT(q3, q1))
circuit.append(cirq.CZ(q0, q1))
circuit.append(cirq.H(q0))
circuit.append(cirq.CNOT(q3, q2))
circuit.append(cirq.ry(rads=3.14/8)(q2))
circuit.append(cirq.ry(rads=-3.14/8)(q3))
circuit.append(cirq.CNOT(q0, q3))
circuit.append(cirq.CNOT(q0, q2))
circuit.append(cirq.ry(rads=-3.14/8)(q3))
circuit.append(cirq.ry(rads=3.14/8)(q2))
circuit.append(cirq.CNOT(q1, q2))
circuit.append(cirq.CNOT(q2, q3))
circuit.append(cirq.ry(rads=-3.14/8)(q2))
circuit.append(cirq.ry(rads=3.14/8)(q3))
circuit.append(cirq.CNOT(q0, q3))
circuit.append(cirq.CZ(q0, q3))
circuit.append(cirq.ry(rads=-3.14/8)(q2))
circuit.append(cirq.ry(rads=3.14/8)(q3))
circuit.append(cirq.CNOT(q3, q2))
circuit.append(cirq.CNOT(q1, q3))
circuit.append(cirq.H(q3))
circuit.append(cirq.CNOT(q3, q1))
circuit.append(cirq.CNOT(q1, q0))

print(circuit)


# In[ ]:





# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange


# In[2]:


def f(t):
    return t ** 10 - 1 * t ** 9

ts = np.linspace(0, 1, 6)
fs = f(ts)
print(ts)


# In[3]:


p4 = lagrange(ts[1:-1], fs[1:-1])
p51 = lagrange(ts[:-1], fs[:-1])
p52 = lagrange(ts[1:], fs[1:])
p6 = lagrange(ts, fs)
ts_cont = np.linspace(0, 1, 100)


# In[4]:


plt.plot(ts_cont, f(ts_cont), ts_cont, p6(ts_cont))


# In[5]:


class LagrangeBary:
    def __init__(self, ts, fs):
        self.ts = np.copy(ts)
        self.fs = np.copy(fs)
        self.compute_ws()
        
    def compute_ws_j(self, j):
        wsj = 1
        for k in range(len(self.ts)):
            if k == j:
                continue
            wsj /= self.ts[j] - self.ts[k]
        return wsj
    
    def compute_ws(self):
        n = len(self.ts)
        ws = np.ones(n)
        for j in range(n):
            ws[j] = self.compute_ws_j(j)
        self.ws = ws
        
    def remove_left(self):
        t = self.ts[0]
        for i in range(len(self.ws)):
            self.ws[i] *= self.ts[i] - t
        self.ts = np.delete(self.ts, 0)
        self.fs = np.delete(self.fs, 0)
        self.ws = np.delete(self.ws, 0)
        
    def add_right(self, ti, fi):
        self.ts = np.append(self.ts, ti)
        self.fs = np.append(self.fs, fi)
        
        for i in range(len(self.ws)):
            self.ws[i] /= self.ts[i] - ti
            
        w_new = self.compute_ws_j(len(self.ts) - 1)
        self.ws = np.append(self.ws, w_new)
        
    def __call__(self, t):
        for i in range(len(self.ts)):
            if self.ts[i] == t:
                return self.fs[i]
        n = len(self.ts)
        nom = sum(self.ws[j] * self.fs[j] / (t - self.ts[j]) for j in range(n))
        denom = sum(self.ws[j] / (t - self.ts[j]) for j in range(n))
        return nom / denom


# In[6]:


p_bary = LagrangeBary(ts[:-1], fs[:-1])


# In[7]:


plt.plot(ts_cont, p51(ts_cont), ts_cont, [p_bary(t) for t in ts_cont])


# In[8]:


p_bary.remove_left()
plt.plot(ts_cont, p4(ts_cont), ts_cont, [p_bary(t) for t in ts_cont])


# In[9]:


p_bary.add_right(ts[-1], f(ts[-1]))
plt.plot(ts_cont, p52(ts_cont), ts_cont, [p_bary(t) for t in ts_cont])


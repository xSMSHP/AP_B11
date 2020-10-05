import re
import numpy as np
import matplotlib.pyplot as plt

class Pair:
    def __init__(self, U):
        self.U = U
        self.short = False
        self.r11a = []
        self.r11i = []
        self.r12a = []
        self.r12i = []
        self.r21 = []
        self.r22 = []
    
    def average(self):
        av = []
        av.append(0)
        for r in self.r11a:
            av[0] += r
        av[0] /= len(self.r11a)
        av.append(0)
        for r in self.r11i:
            av[1] += r
        av[1] /= len(self.r11i)
        av.append(0)
        for r in self.r12a:
            av[2] += r
        av[2] /= len(self.r12a)
        av.append(0)
        for r in self.r12i:
            av[3] += r
        av[3] /= len(self.r12i)
        if len(self.r21) * len(self.r22) != 0:
            av.append(0)
            for r in self.r21:
                av[4] += r
            av[4] /= len(self.r21)
            av.append(0)
            for r in self.r22:
                av[5] += r
            av[5] /= len(self.r22)
        return av
    
    def variance(self):
        var = []
        av = self.average()
        var.append(0)
        for r in self.r11a:
            var[0] += (av[0] - r) ** 2
        var[0] = np.sqrt(var[0] / (len(self.r11a) - 1))
        var.append(0)
        for r in self.r11i:
            var[1] += (av[1] - r) ** 2
        var[1] = np.sqrt(var[1] / (len(self.r11i) - 1))
        var.append(0)
        for r in self.r12a:
            var[2] += (av[2] - r) ** 2
        var[2] = np.sqrt(var[2] / (len(self.r12a) - 1))
        var.append(0)
        for r in self.r12i:
            var[3] += (av[3] - r) ** 2
        var[3] = np.sqrt(var[3] / (len(self.r12i) - 1))
        if len(self.r21) * len(self.r22) != 0:
            var.append(0)
            for r in self.r21:
                var[4] += (av[4] - r) ** 2
            var[4] = np.sqrt(var[4] / (len(self.r21) - 1))
            var.append(0)
            for r in self.r22:
                var[5] += (av[5] - r) ** 2
            var[5] = np.sqrt(var[5] / (len(self.r22) - 1))
        return var
    
    def output(self):
        av = self.average()
        out = []
        out.append(np.sqrt((av[0] ** 2 + av[1] ** 2) / 2))
        out.append(np.sqrt((av[2] ** 2 + av[3] ** 2) / 2))
        if self.short:
            out.append(av[4])
            out.append(av[5])
        return out

class A1:
    def __init__(self, file, name):
        self.name = name
        self.P = []
        
        for line in file:
            l = re.split(',', line)
            ind = -1
            for i, p in enumerate(self.P):
                if p.U == float(l[0]):
                    ind = i
            if ind == -1:
                self.P.append(Pair(float(l[0])))
            self.add(l, ind)
                
    def add(self, l, i):
        self.P[i].r11a.append(float(l[1]))
        self.P[i].r11a.append(float(l[3]))
        self.P[i].r11i.append(float(l[2]))
        self.P[i].r11i.append(float(l[4]))
        self.P[i].r12a.append(float(l[5]))
        self.P[i].r12a.append(float(l[7]))
        self.P[i].r12i.append(float(l[6]))
        self.P[i].r12i.append(float(l[8]))
        if int(l[9]) != -1:
            if not self.P[i].short:
                self.P[i].short = True
            self.P[i].r21.append(float(l[9]))
            self.P[i].r21.append(float(l[10]))
            self.P[i].r22.append(float(l[11]))
            self.P[i].r22.append(float(l[12]))
    
    def variance(self):
        var = 0
        count = 0
        for p in self.P:
            p_var = p.variance()
            for v in p_var:
                var += v
                count += 1
        return var / count
    
    def plot(self):
        r11 = []
        r12 = []
        r21 = []
        r22 = []
        λ = []
        λ_short = []
        for p in self.P:
            λ.append(1225 / np.sqrt(p.U * 10 ** 3))
            if p.short:
                λ_short.append(1225 / np.sqrt(p.U * 10 ** 3))
            o = p.output()
            r11.append(o[0])
            r12.append(o[1])
            if p.short:
                r21.append(o[2])
                r22.append(o[3])
        
        plt.figure()
        
        coef1 = np.polyfit(λ, r11, 1)
        poly1d_fn1 = np.poly1d(coef1) 
        coef2 = np.polyfit(λ, r12, 1)
        poly1d_fn2 = np.poly1d(coef2) 
        print(poly1d_fn1, poly1d_fn2)
        print(127 * 10 ** (-3) / (coef1[0] * 10 ** 9))
        print(127 * 10 ** (-3) / (coef2[0] * 10 ** 9))
        
        plt.errorbar(λ, r11, yerr = self.variance(), marker = 's', color = 'red', linestyle = '')
        plt.errorbar(λ, r12, yerr = self.variance(), marker = 'd', color = 'blue', linestyle = '')
        plt.plot(np.linspace(12, 16, 100), poly1d_fn1(np.linspace(12, 16, 100)), 'r-')
        plt.plot(np.linspace(12, 16, 100), poly1d_fn2(np.linspace(12, 16, 100)), 'b-')
        plt.xlabel(r'$\lambda$ [pm]')
        plt.ylabel('r [mm]')
        plt.legend([r'Ausgleichsgerade $d_1$', r'Ausgleichsgerade $d_2$', r'Messwerte $d_1$', r'Messwerte $d_2$'])
        plt.grid(True)
        plt.autoscale(True)
        plt.savefig(self.name + '1.png', dpi = 900)
        plt.show()
        
        coef3 = np.polyfit(λ_short, r21, 1)
        poly1d_fn3 = np.poly1d(coef3) 
        coef4 = np.polyfit(λ_short, r22, 1)
        poly1d_fn4 = np.poly1d(coef4) 
        print(poly1d_fn3, poly1d_fn4)
        print(127 * 10 ** (-3) * 2 / (coef3[0] * 10 ** 9))
        print(127 * 10 ** (-3) * 2 / (coef4[0] * 10 ** 9))
        
        plt.figure()
        plt.errorbar(λ_short, r21, yerr = self.variance(), marker = 'd', color = 'red', linestyle = '')
        plt.errorbar(λ_short, r22, yerr = self.variance(), marker = 's', color = 'blue', linestyle = '')
        plt.plot(np.linspace(12, 14, 100), poly1d_fn3(np.linspace(12, 14, 100)), 'r-')
        plt.plot(np.linspace(12, 14, 100), poly1d_fn4(np.linspace(12, 14, 100)), 'b-')
        plt.xlabel(r'$\lambda$ [pm]')
        plt.ylabel('r [mm]')
        plt.legend([r'Ausgleichsgerade $d_2$', r'Ausgleichsgerade $d_1$', r'Messwerte $d_2$', r'Messwerte $d_1$'])
        plt.grid(True)
        plt.autoscale(True)
        plt.savefig(self.name + '2.png', dpi = 900)
        plt.show()

a1 = A1(open('b11.txt', 'r'), 'b11')
a1.plot()
print('Δr = ' + str(a1.variance()))

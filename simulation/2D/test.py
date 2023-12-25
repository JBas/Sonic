import numpy as np
from skimage.io import imsave, imshow, show
import skimage.draw
import matplotlib.pyplot as plt
from scipy.special import j1

c = 343*1000
f = 40000
L = c/f
k = 2*np.pi/L

P0 = 1
A = 15
a = 4
diam = int(2*a)
N = 15

width = 256
width_2 = width//2
height = 256

def unit_vector(vector):
    return vector / np.linalg.norm(vector)

def calcTheta(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)

    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

class Transducer:
    def __init__(self, pt, isTop=False):
        self.pt = pt
        self.phi = 0
        self.isTop = isTop
        self.norm =  np.array([0, 1]) if isTop else np.array([0, -1])
        pass

    def calcPhase(self, pt):
        delta = pt - self.pt
        d = np.linalg.norm(delta)
        phase = -d*k
        if (self.isTop):
            phase += np.pi
        self.phi = np.mod(phase, 2*np.pi)

    def Df(theta):
        # return 2*j1(k*a*np.sin(theta)) / k*a*np.sin(theta)
        return np.sinc(k*a*np.sin(theta))

    def P(self, pt):
        delta = pt-self.pt
        d = np.linalg.norm(delta)
        if d == 0:
            return 0, 0
        theta = calcTheta(pt+self.pt, self.norm+self.pt)
        real = P0*A*Transducer.Df(theta)*np.cos(self.phi + k*d)
        im = P0*A*Transducer.Df(theta)*np.sin(self.phi + k*d)
        return real, im
    

if __name__=="__main__":

    frame = 0

    transducers = []
    for i in range(-N//2 + 1, N//2 + 1):
        transducers.append(Transducer(np.array([width_2 + i*diam, 0]), True))

    for i in range(-N//2 + 1, N//2 + 1):
        transducers.append(Transducer(np.array([width_2 + i*diam, height-1]), False))

    for _x in range(0, width, 16):
        img = np.zeros([height, width, 3])
        focus = np.array([_x, height//2])

        for transducer in transducers:
            transducer.calcPhase(focus)

        for y in range(1, height-1):
            for x in range(width):
                real = 0
                im = 0
                for transducer in transducers:
                    _real, _im = transducer.P(np.array([x, y]))
                    real += _real
                    im += _im
                power = np.sqrt(real*real + im*im)
                img[y, x] = [power, power, power]


        img = img-np.min(img)
        img = 255*img/np.max(img)

        rr, cc = skimage.draw.disk(tuple([focus[1], focus[0]]), 4.5, shape=img.shape)
        img[rr, cc] = [255, 0, 0]

        for transducer in transducers:
            x = transducer.pt[0]
            y = transducer.pt[1]
            rr, cc = skimage.draw.disk(tuple([y, x]), a, shape=img.shape)
            img[rr, cc] = [0, 255, 0]

            if transducer.isTop:
                rr, cc = skimage.draw.line(y, x, y + 25, x)
            else:
                rr, cc = skimage.draw.line(y, x, y - 25, x)
            img[rr, cc] = [0, 0, 255]

        img = img.astype(np.uint8)

        imsave(f"frames/aframe_{frame}.png", img)
        frame += 1
        # imshow(img)
        # show()
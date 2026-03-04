import numpy as np
from numpy import linalg as LA

#TODO
'''
Decide if delta quat from gyro is a class method or its own function
Implement rotate and to_euler/to_rotation matrix

'''

class quaternion(np.ndarray):
    def __new__(cls, w, x, y, z):
        
        arr = np.array([w, x, y, z], dtype=np.float64)
        
        return arr.view(cls)
    
    def __array_finalize__(self, obj):
        if obj is None: return
    
    @property
    def w(self):
        return self[0]
    @property
    def x(self):
        return self[1]
    @property
    def y(self):
        return self[2]
    @property
    def z(self):
        return self[3]
    
    @w.setter
    def w(self, value):
        self[0] = value

    @x.setter
    def x(self, value):
        self[1] = value
    
    @y.setter
    def y(self, value):
        self[2] = value
    
    @z.setter   
    def z(self, value):
        self[3] = value

    
    def norm(self):
       return LA.norm(self)
    
    def unit_quat(self):
        return self / LA.norm(self)
    
    def conjugate(self): 
        x = self.x * -1
        y = self.y * -1
        z = self.z * -1
        quat = quaternion(self.w, x, y, z)
        return quat
    
    def __mul__(self, other):
        w = (self.w*other.w) - (self.x*other.x) - (self.y*other.y) - (self.z*other.z)
        x = (self.w*other.x) + (self.x*other.w) + (self.y*other.z) - (self.z*other.y)
        y = (self.w*other.y) - (self.x*other.z) + (self.y*other.w) + (self.z*other.x)
        z = (self.w*other.z) + (self.x*other.y) - (self.y*other.x) + (self.z*other.w)

        return quaternion(w, x, y, z)
    
    def rotate(self):
        raise NotImplementedError
    

    
    def deltaQuatFromGyro(self, dt):
        epsilon = 1e-12 #some aribtrary value to avoid floating point error
    
        theta = self.norm() * dt

        if theta < epsilon:
            return quaternion(1, 0, 0, 0)

        u = self.unit_quat()

        s = np.sin(theta*0.5)
        c = np.cos(theta*0.5)

        dq = quaternion(c, u.x*s, u.y*s, u.z*s)

        dq = dq.unit_quat()

        return dq
        




   


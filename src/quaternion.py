import torch

class TensorMixin:
    def norm(self):
        return torch.linalg.vector_norm(self._data)
    def unit(self):
        return type(self)(self._data / self.norm())
    def __repr__(self):
        return f"{type(self)}: {self._data}"
    

class GyroTensor(TensorMixin):
    def __init__(self, data: torch.tensor):
        assert data.shape[-1] == 3, 'Gyro Dim must be 3'
        self._data = torch.as_tensor(data , dtype=torch.float64)
    
    @property
    def x(self):
        return self._data[0]
    
    @x.setter
    def x(self, value):
        self._data[0] = value

    @property
    def y(self):
        return self._data[1]
    
    @y.setter
    def y(self, value):
        self._data[1] = value

    @property
    def z(self):
        return self._data[2]
    
    @z.setter
    def z(self, value):
        self._data[2] = value

    def deltaQuatFromGyro(self, dt):
        epsilon = 1e-12 #some aribtrary value to avoid floating point error
    
        theta = self.norm() * dt

        if theta < epsilon:
            return QuaternionTensor(torch.tensor([1,0,0,0]))

        u = self.unit()

        s = torch.sin(theta*0.5)
        c = torch.cos(theta*0.5)

        dq = QuaternionTensor(torch.tensor([c, u.x*s, u.y*s, u.z*s]))

        dq = dq.unit()

        return dq


class QuaternionTensor(TensorMixin):
    def __init__(self, data: torch.tensor):
        assert data.shape[-1] == 4, 'Quaternion Dim must be 4'
        self._data = torch.as_tensor(data, dtype=torch.float64)
    
    @property
    def w(self):
        return self._data[0]
    
    @w.setter
    def w(self, value):
        self._data[0] = value

    @property
    def x(self):
        return self._data[1]
    
    @x.setter
    def x(self, value):
        self._data[1] = value

    @property
    def y(self):
        return self._data[2]
    
    @y.setter
    def y(self, value):
        self._data[2] = value

    @property
    def z(self):
        return self._data[3]
    
    @z.setter
    def z(self, value):
        self._data[3] = value

    def conjugate(self):
        x = self.x * -1
        y = self.y * -1
        z = self.z * -1
        
        quat = QuaternionTensor(torch.stack([self.w, x, y, z]))
        return quat 

    def rotate(self):
        raise NotImplementedError

    def __mul__(self, other):
        w = (self.w*other.w) - (self.x*other.x) - (self.y*other.y) - (self.z*other.z)
        x = (self.w*other.x) + (self.x*other.w) + (self.y*other.z) - (self.z*other.y)
        y = (self.w*other.y) - (self.x*other.z) + (self.y*other.w) + (self.z*other.x)
        z = (self.w*other.z) + (self.x*other.y) - (self.y*other.x) + (self.z*other.w)

        return QuaternionTensor(torch.stack([w,x,y,z]))
    
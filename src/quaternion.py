import torch


class QuaternionTensor:
    def __init__(self, data: torch.tensor):
        assert data.shape[-1] == 4, 'Dim must be 4'

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

    def norm(self):
        return torch.linalg.vector_norm(self._data)

    def unit_quat(self):
        return QuaternionTensor(self._data / self.norm())

    def conjugate(self):
        x = self.x * -1
        y = self.y * -1
        z = self.z * -1

        tensor = torch.stack([self.w, x, y, z])
        quat = QuaternionTensor(tensor)
        return quat 

    def rotate(self):
        raise NotImplementedError

    def __mul__(self, other):
        w = (self.w*other.w) - (self.x*other.x) - (self.y*other.y) - (self.z*other.z)
        x = (self.w*other.x) + (self.x*other.w) + (self.y*other.z) - (self.z*other.y)
        y = (self.w*other.y) - (self.x*other.z) + (self.y*other.w) + (self.z*other.x)
        z = (self.w*other.z) + (self.x*other.y) - (self.y*other.x) + (self.z*other.w)

        tensor = torch.stack([w,x,y,z])
        return QuaternionTensor(tensor)
    
    def __repr__(self):
        return f"QuaternionTensor: {self._data}"
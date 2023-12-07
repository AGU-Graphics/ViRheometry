import plotly.express as px
import numpy as np

class ConstMeta(type):
    def __setattr__(self, name, value):
        if name in self.__dict__:
            raise TypeError(f'Can\'t rebind const ({name})')
        else:
            self.__setattr__(name, value)

   
class CGS(metaclass=ConstMeta):
    COLORS = px.colors.qualitative.Plotly
    MIN_LOSS_THR = 0.0
    MAX_LOSS_THR = 1000.0


    RESOLUTION_ETA = 30
    RESOLUTION_N = 30
    RESOLUTION_SIGMA_Y = 30

    RESOLUTION_X = 400
    RESOLUTION_Y = 400

    MIN_ETA = 0.001
    MAX_ETA = 300.0

    MIN_N = 0.3
    MAX_N = 1.0

    MIN_SIGMA_Y = 0.0
    MAX_SIGMA_Y = 400.0

    MIN_W = 2.0
    MAX_W = 7.0

    MIN_H = 2.0
    MAX_H = 7.0
        
    MIN_Wmm = 20
    MAX_Wmm = 70

    MIN_Hmm = 20
    MAX_Hmm = 70

    GRAVITY = 981.0

    WEIGHTS_STEPS = 100

    extent_eta    = np.array([MIN_ETA, MAX_ETA])
    extent_n      = np.array([MIN_N, MAX_N])
    extent_sigmaY = np.array([MIN_SIGMA_Y, MAX_SIGMA_Y])
    mmToCurrentUnit = 1e-1


class MKS(metaclass=ConstMeta):
    COLORS = px.colors.qualitative.Plotly
    MIN_LOSS_THR = CGS.MIN_LOSS_THR
    MAX_LOSS_THR = CGS.MAX_LOSS_THR

    RESOLUTION_ETA = CGS.RESOLUTION_ETA
    RESOLUTION_N = CGS.RESOLUTION_N
    RESOLUTION_SIGMA_Y = CGS.RESOLUTION_SIGMA_Y

    RESOLUTION_X = CGS.RESOLUTION_X
    RESOLUTION_Y = CGS.RESOLUTION_Y

    MIN_ETA = CGS.MIN_ETA * 1e-1
    MAX_ETA = CGS.MAX_ETA * 1e-1

    MIN_N = CGS.MIN_N
    MAX_N = CGS.MAX_N

    MIN_SIGMA_Y = CGS.MIN_SIGMA_Y * 1e-1
    MAX_SIGMA_Y = CGS.MAX_SIGMA_Y * 1e-1

    MIN_W = CGS.MIN_W * 1e-2
    MAX_W = CGS.MAX_W * 1e-2
    MIN_H = CGS.MIN_H * 1e-2
    MAX_H = CGS.MAX_H * 1e-2

    MIN_Wmm = CGS.MIN_Wmm
    MAX_Wmm = CGS.MAX_Wmm
    MIN_Hmm = CGS.MIN_Hmm
    MAX_Hmm = CGS.MAX_Hmm

    GRAVITY = CGS.GRAVITY * 1e-2
    mmToCurrentUnit = 1e-3
    WEIGHTS_STEPS = CGS.WEIGHTS_STEPS

    extent_eta    = np.array([MIN_ETA, MAX_ETA])
    extent_n      = np.array([MIN_N, MAX_N])
    extent_sigmaY = np.array([MIN_SIGMA_Y, MAX_SIGMA_Y])
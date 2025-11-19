import Parameters as param
import IsentropicTools as IT
import numpy as np

T0 = param.T_combustion
g = param.g
Rs = param.Rs


def LocalTemperature(wall_x, wall_y):
    T_list = []
    for i in range(len(wall_x)):
        if wall_x[i] < 0:
                flow_regime = 'subsonic'
        else:
            flow_regime = 'supersonic'
        
        AR = wall_y[i]**2 / np.min(wall_y)**2
        M_local = IT.AreaRatioInverse(AR, g, flow_regime)
        T_local = IT.Temperature(T0, g, M_local)
        T_list.append(T_local)
    nan_indices = np.where(np.isnan(T_list))[0]
    for j in nan_indices:
         T_list[j] = T_list[j - 1] / 2 + T_list[j + 1] / 2
              
    return np.array(T_list)
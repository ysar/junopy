import numpy as np
import datetime as dt
import pandas as pd



def _construct_F_anode():
    """
    Constructs the constant table of F(A) where F corresponds to the equations 
    in 
    and A refers to anode.
    """
    df = pd.DataFrame(columns=['c0', 'c1', 'c2', 'b0', 'b1', 'd', 'f'])

    df['c0'] = np.full((32,), 3.12e-5)
    df['c1'] = np.full((32,), -5.3e-8)
    df['c2'] = np.full((32,), -6.92e-9)
    df['b0'] = np.concatenate((np.full((16,), 1.088), np.full((16,), 1.048)))
    df['b1'] = np.concatenate((
        np.full((16,), -9.766e-3), np.full((16,), -5.372e-3)))
    df['f'] = np.full((32,), 1.)
    df['d'] = np.concatenate((np.full((16,), 16), np.full((16,), 32)))
    df.loc[[0, 15, 16, 31], 'f'] = 0.745
    
    return df


def _calc_deflection_effect(delta):
    """
    Assume that delta has shape (Anode, Energy, time). 
    This is NOT supposed to be used for LRS data.
    """
    assert (delta.shape[0] == 32)
    dfConsts = _construct_F_anode()

    anode = dfConsts.index

    D_delta = delta * dfConsts['c1'] / dfConsts['c0']\
            + delta * dfConsts['c2'] / dfConsts['c0']
    
    F_A = dfConsts['f'] * (
        dfConsts['b0'] + dfConsts['b1'] * (dfConsts['d'] - anode.T))
        
    return (F_A * D_delta).T


def get_Geff_electrons(energies, times, elevations, sensor='BOTH'):

    rowindex = ['060', '180']
    columns = ['tau0', 'tau1', 'Emax', 'rho', 'n', 'k', 'epsilonmax']
    paramCases = [
        pd.DataFrame(index=rowindex, columns=columns).fillna(0.), 
        pd.DataFrame(index=rowindex, columns=columns).fillna(0.), 
        pd.DataFrame(index=rowindex, columns=columns).fillna(0.)]
    
    # Case 0 : Before any corrections
    paramCases[0]['tau0'] = [1.04, 1.16]
    paramCases[0]['Emax'] = [1811, 1837]
    paramCases[0]['rho'] = [2.22, 3.104]
    paramCases[0]['n'] = 1.67
    paramCases[0]['k'] = [4.033, 3.352]
    paramCases[0]['epsilonmax'] = [0.923, 1.137]

    # Case 1 : 2016-07-04 to 2019-05-21
    paramCases[1]['tau0'] = [1.04, 1.16]
    paramCases[1]['tau1'] = [-0.062, -0.079]
    paramCases[1]['Emax'] = [870, 1140]
    paramCases[1]['rho'] = [1.33, 2.67]
    paramCases[1]['n'] = 1.67
    paramCases[1]['k'] = [4.0, 2.86]
    paramCases[1]['epsilonmax'] = [0.95, 1.21]

    # Case 2 : 2019-05-22 to -
    paramCases[2]['tau0'] = [1.04, 1.16]
    paramCases[2]['tau1'] = [-0.062, -0.079]
    paramCases[2]['Emax'] = [870, 1140]
    paramCases[2]['rho'] = [1.33, 2.67]
    paramCases[2]['n'] = 1.67
    paramCases[2]['k'] = [4.0, 2.86]
    paramCases[2]['epsilonmax'] = [0.95, 1.21]

    
    # shape(nDir, nE, nTime)
    G_simulation = 0.8 * 0.9 * (2.37e-4 + 9.17e-11 * energies)
    
    G_eff = G_simulation
    
    # Consider anode dependence if in HRS
    if energies.shape[0] == 32: 
        G_eff *= _calc_deflection_effect(elevations)

    delta = np.array([
        (dt.datetime(2016, 7, 4) - x) / dt.timedelta(year=1) \
            if x > dt.datetime(2016, 7, 4) else 0. for x in times ])
    
    print(delta)

    effCase = np.zeros(times.shape)
    for i, t in enumerate(times):
        if t < dt.datetime(2016, 7, 3):
            effCase[:, :, i] = 0
        elif t < dt.datetime(2019, 5, 22):
            effCase[:, :, i] = 1
        else:
            effCase[:, :, i] = 2

    tau = np.zeros(energies)
    eff = np.zeros(energies)
    for i in [0, 1, 2]:
        iCase = (effCase == 0).nonzero()

        if energies.shape[0] == 32:    # Only two sensors

            col_060 = paramCases[iCase].loc['060', :]
            col_180 = paramCases[iCase].loc['180', :]

            tau[:16, :, iCase] = \
                col_060['tau0'] + col_060['tau1'] * delta[iCase]
            
            tau[16:, :, iCase] = \
                col_180['tau0'] + col_180['tau1'] * delta[iCase]
            
            eff[:16 ,:, iCase] = \
                col_060['epsilonmax'] \
                * (1 - np.exp(
                    -col_060['k'] 
                    * (energies[:16, :, iCase] / col_060['Emax'])**(1 - col_060['n'])
                    * (1 - np.exp(-col_060['rho'] 
                                  * (energies[:16, :, iCase])**col_060['n'])
                        )
                    )
                )
            
            eff[16: ,:, iCase] = \
                col_180['epsilonmax'] \
                * (1 - np.exp(
                    -col_180['k'] 
                    * (energies[16:, :, iCase] / col_180['Emax'])**(1 - col_180['n'])
                    * (1 - np.exp(-col_180['rho'] 
                                  * (energies[16:, :, iCase])**col_180['n'])
                        )
                    )
                )
                
            
        else:
            
            col = paramCases[iCase].loc[sensor, :]

            tau[:, :, iCase] = col['tau0'] + col['tau1'] * delta[iCase]
            
            eff[: ,:, iCase] = \
                col['epsilonmax'] \
                * (1 - np.exp(
                    -col['k'] 
                    * (energies[:, :, iCase] / col['Emax'])**(1 - col['n'])
                    * (1 - np.exp(-col['rho'] 
                                  * (energies[:, :, iCase])**col['n'])
                        )
                    )
                )
            
    Geff *= tau * eff

    return Geff




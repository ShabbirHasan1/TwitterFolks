# https://gist.github.com/jace48/34cac8e8b5275f7a33835fecd96821d5
# By `jace48`

import math
from scipy.stats import norm


def BSM_withAdjoints(S0, r, y, sig, K, T):
    #Evaluation
    
    sqrtT = math.sqrt(T)
    
    df = math.exp(-r * T)
    F = S0 * math.exp((r - y) * T)
    std = sig * sqrtT
    d = math.log(F / K) / std
    d1, d2 = d + 0.5 * std, d - 0.5 * std

    nd1, nd2 = norm.cdf(d1), norm.cdf(d2)
    Call_P = df * (F * nd1 - K * nd2)
    Put_P = df * (K * (1 - nd2) - F * (1 - nd1))
    
    # Adjoint calculation
    
    v_ = 1.0
    
    df_ = v_ * (F * nd1 - K * nd2)
    F_ = v_ * df * nd1
    nd1_ = v_ * df * F
    if K:
        K_ = - v_ * df * nd2
    nd2_ = - v_ * df * K
    
    d2_ = nd2_ * norm.pdf(d2)
    d1_ = nd1_ * norm.pdf(d1)
    
    d_ = d2_
    std_ = - 0.5 * d2_
    d_ += d1_
    std_ += 0.5 * d1_
    
    F_ += d_ / (F * std)
    K_ -= d_ / (K * std)
    std_ -= d_ * d / std
    
    sig_ = std_ * sqrtT
    T_ = 0.5 * std_ * sig / sqrtT
    
    S0_ = F_ * F / S0
    r_ = F_ * T * F if r == 0 else 0 
    
    y_ = - F_ * T * F if y == 0 else 0
    T_ += F_ * (r - y) * F
    

    r_ += - df_ * df * T  if r == 0 else 0
    T_ += - df_ * df * r
    
    return (Call_P,Put_P), (S0_, r_, y_, sig_, K_, T_)



if __name__ == "__main__":
    #Example Fwd=22045, Rate & Yeild=0,  Vol=0.1352, Strike=22050, DTE=6 (0.016438356)
    OptionPrice,adjoints= BSM_withAdjoints (22045, 0, 0, 0.1352, 22050, 0.016438356)    
    
    print("Call Option Price:", OptionPrice[0])
    print("Put Option Price:", OptionPrice[1])

    print("S0/Delta:", adjoints[0])
    print("r/Rho:", adjoints[1])
    print("y/Dividend yeild:", adjoints[2])
    print("sig/Vega:", (adjoints[3]))
    print("K/Strike:", adjoints[4])
    print("T/Theta:", (adjoints[5]))  
    

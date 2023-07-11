
import numpy as np
def MatRestrULA(M,J,fs,NumAng):
    
    d = 0.04287
    Ts = 1/fs
    c = 343.0
    mu = d/(c*Ts)
    retardo = 25*Ts
    fmin = 100.0
    fmax = 4000.0
    C_list = []
    Bk_list = []
    fr_list = []
    Xc = np.zeros(M*J)
    Xs = np.zeros(M*J)
    m = np.array(range(M))

    ang_list = np.linspace(-np.pi/2,np.pi/2,NumAng,endpoint=True)

    for a in ang_list:

        r = int(np.ceil(2*(fmax-fmin) * ((M-1)*mu*np.abs(np.sin(a))+(J-1))*Ts + 1))+5
        omega = np.linspace(2*np.pi*fmin*Ts, 2*np.pi*fmax*Ts, 10*r)
        
        l = len(omega)
        C = np.zeros((M*J,2*l),float) 
        f = np.ones((2*l),float)       

        idx = -1
        for w in omega:
            for j in range(J):
                Xc[j*M:(j+1)*M] = np.cos(w*(m*mu*np.sin(a)+j))
                Xs[j*M:(j+1)*M] = np.sin(w*(m*mu*np.sin(a)+j))
                pass
            Xs[0] = 1
            
            idx += 1
            C[:,idx] = Xc
            f[idx] = np.cos((w/Ts)*retardo)     
            C[:,idx+l] = Xs
            f[idx+l] = np.sin((w/Ts)*retardo)
            pass

        U, s, Vh = np.linalg.svd(C)

        C_list.append(U[:,:r])
        Bk_list.append(U[:,r:])
        Sr = np.diag(s[:r])
        Vhr = Vh[:r,:]
        fr = np.dot(np.linalg.inv(Sr),np.dot(Vhr,f))
        fr_list.append(fr)
        pass
    return C_list,fr_list,Bk_list

def MatRestrUCA(M,J,fs,NumAng):
    
    Ts = 1/fs
    c = 343.0
    fmin = 100.0
    fmax = 4000.0
    lamb = c/fmax
    rc = lamb/(4*np.sin(np.pi/M))
    mu = rc/(c*Ts)
    retardo = 25*Ts
    C_list = []
    Bk_list = []
    fr_list = []
    Xc = np.zeros(M*J)
    Xs = np.zeros(M*J)

    ang_list = np.linspace(0,2*np.pi,NumAng,endpoint=False)
    sensor_angs = np.linspace(0,2*np.pi, M,endpoint=False)

    for a in ang_list:

        r = int(np.ceil(2*(fmax-fmin) * ((M-1)*mu*np.abs(np.sin(a))+(J-1))*Ts + 1))+5
        omega = np.linspace(2*np.pi*fmin*Ts, 2*np.pi*fmax*Ts, 10*r)
        
        l = len(omega)
        C = np.zeros((M*J,2*l),float)
        f = np.ones((2*l),float) 

        idx = -1
        for w in omega:
            for j in range(J):
                Xc[j*M:(j+1)*M] = np.cos(w*(mu*np.cos(a - sensor_angs)+j))
                Xs[j*M:(j+1)*M] = np.sin(w*(mu*np.cos(a - sensor_angs)+j))
            
            idx += 1
            C[:,idx] = Xc
            f[idx] = np.cos((w/Ts)*retardo)     
            C[:,idx+l] = Xs
            f[idx+l] = np.sin((w/Ts)*retardo)
            pass

        U, s, Vh = np.linalg.svd(C)

        C_list.append(U[:,:r])
        Bk_list.append(U[:,r:])
        Sr = np.diag(s[:r])
        Vhr = Vh[:r,:]
        fr = np.dot(np.linalg.inv(Sr),np.dot(Vhr,f))
        fr_list.append(fr)
        pass
    return C_list,fr_list,Bk_list


def MultiplesDOAs(array, NumAng, N=1, array_mic="ULA"):
    list_indx = []
    doas = []
    indices =[]
    for i in array:
        indx_mx = i.argmax()
        if N>1:
            a_indx = np.where(i >= (i[indx_mx]/10))
            for a in a_indx[0]:
                m1 = (i[a]-i[a-1])/(a-(a-1))
                if a < len(i)-1:
                    m2 = (i[a+1]-i[a])/((a+1)-a)
                else:
                    m2 = m1
                if np.logical_and(m1>0,m2<0):
                    list_indx.append(a)
                    pass
                pass
            vals=i[list_indx][i[list_indx].argsort()[::-1][:N]]
            ind = [np.where(i==elto)[0].tolist()[0] for elto in vals]
            indices.append(ind)
            list_indx = []
            if array_mic == "ULA":
                doas.append(list(map(lambda x: (np.pi/(NumAng)*(x)-np.pi/2)*180/np.pi, ind)))
            elif array_mic == "UCA":
                doas.append(list(map(lambda x: (2*np.pi/(NumAng)*(x))*180/np.pi, ind)))
        else:
            if array_mic == "ULA":
                doas += [(np.pi/(NumAng)*(indx_mx)-np.pi/2)*180/np.pi]
            elif array_mic == "UCA":
                doas += [(2*np.pi/(NumAng)*(indx_mx))*180/np.pi]
            indices += [indx_mx]
    return (doas, indices)


def getIdx(doas, NumAng, array_mic = 'ULA'):
    index = []
    if array_mic == 'ULA':
        for i in doas:
            index.append(int(np.round(i+90)))
    elif array_mic == 'UCA':
        for i in doas:
            index.append(int(np.round(((360+i)-360*0**np.heaviside(-i,0))/2)))
    return index
    
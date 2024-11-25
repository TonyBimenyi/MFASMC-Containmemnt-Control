import numpy as np
import matplotlib.pyplot as plt

#Defining Parameters

beta = -2
rho = 1
eta = 40
lambda_ = 10
mu = 1
epsilon = 10**-5
alpha1 = 1
alpha2 = 1.5
T = 0.1
gamma1 = 0.15
gamma2 = 0.15
gamma3 = 0.45
gamma4 = 0.45
rT = 1024
m = 350
L = 500

# Generate constant desired trajectory
yd = 0.5 * (1 + np.sign(np.sin(np.linspace(0, 2 * np.pi, L + 1))))  # Square wave

# Initialize arrays
phi1 = np.zeros((L, 1))
phi2 = np.zeros((L, 1))
phi3 = np.zeros((L, 1))
phi4 = np.zeros((L, 1))

mfa1 = np.zeros((L, 1))
mfa2 = np.zeros((L, 1))
mfa3 = np.zeros((L, 1))
mfa4 = np.zeros((L, 1))

sm1 = np.zeros((L, 1))
sm2 = np.zeros((L, 1))
sm3 = np.zeros((L, 1))
sm4 = np.zeros((L, 1))

u1 = np.zeros((L, 1))
u2 = np.zeros((L, 1))
u3 = np.zeros((L, 1))
u4 = np.zeros((L, 1))

y1 = np.zeros((L + 1, 1))
y2 = np.zeros((L + 1, 1))
y3 = np.zeros((L + 1, 1))
y4 = np.zeros((L + 1, 1))

e1 = np.zeros((L + 1, 1))
e2 = np.zeros((L + 1, 1))
e3 = np.zeros((L + 1, 1))
e4 = np.zeros((L + 1, 1))

w5 = np.zeros((L, 1))
w6 = np.zeros((L, 1))

si1 = np.zeros((L, 1))
si2 = np.zeros((L, 1))
si3 = np.zeros((L, 1))
si4 = np.zeros((L, 1))

ss1= np.zeros((L, 1))
ss2= np.zeros((L, 1))
ss3= np.zeros((L, 1))
ss4= np.zeros((L, 1))


# Set w5 and w6 according to conditions
w5[:165] = 1.4
w5[165:330] = 1.6
w5[330:] = 1.3

w6[:165] = 0.7
w6[165:330] = 1.2
w6[330:] = 1.1





# Simulation loop

for k in range (1, L-1):
    if k == 1:
        phi1[1] = 1
        phi2[1] = 1
        phi3[1] = 1
        phi4[1] = 1
    elif k == 2:
        phi1[k] = phi1[k - 1] + eta * u1[k - 1] / (mu + u1[k - 1]**2) * (y1[k] - y1[k - 1] - phi1[k - 1] * u1[k - 1])
        phi2[k] = phi2[k - 1] + eta * u2[k - 1] / (mu + u2[k - 1]**2) * (y2[k] - y2[k - 1] - phi2[k - 1] * u2[k - 1])
        phi3[k] = phi3[k - 1] + eta * u3[k - 1] / (mu + u3[k - 1]**2) * (y3[k] - y3[k - 1] - phi3[k - 1] * u3[k - 1])
        phi4[k] = phi4[k - 1] + eta * u4[k - 1] / (mu + u4[k - 1]**2) * (y4[k] - y4[k - 1] - phi4[k - 1] * u4[k - 1])
    else:
        phi1[k] = phi1[k - 1] + (eta * (u1[k - 1] - u1[k - 2]) / (mu + (abs(u1[k - 1] - u1[k - 2]))**2)) * (y1[k] - y1[k - 1] - phi1[k - 1] * (u1[k - 1] - u1[k - 2]))
        phi2[k] = phi2[k - 1] + (eta * (u2[k - 1] - u2[k - 2]) / (mu + (abs(u2[k - 1] - u2[k - 2]))**2)) * (y2[k] - y2[k - 1] - phi2[k - 1] * (u2[k - 1] - u2[k - 2]))
        phi3[k] = phi3[k - 1] + (eta * (u3[k - 1] - u3[k - 2]) / (mu + (abs(u3[k - 1] - u3[k - 2]))**2)) * (y3[k] - y3[k - 1] - phi3[k - 1] * (u3[k - 1] - u3[k - 2]))
        phi4[k] = phi4[k - 1] + (eta * (u4[k - 1] - u4[k - 2]) / (mu + (abs(u4[k - 1] - u4[k - 2]))**2)) * (y4[k] - y4[k - 1] - phi4[k - 1] * (u4[k - 1] - u4[k - 2]))


    si1[k] = y2[k] - 2 * y1[k] + w5[k]
    si2[k] = y3[k] - y2[k]
    si3[k] = y4[k] - 2 * y3[k] + y1[k]
    si4[k] = y2[k] - 2 * y4[k] + w6[k]



    # ss1[k] = si1[k] + alpha1*si1[k-1] + alpha2*np.sign(si1[k-1])
    # ss2[k] = si2[k] + alpha1*si2[k-1] + alpha2*np.sign(si2[k-1])
    # ss3[k] = si3[k] + alpha1*si3[k-1] + alpha2*np.sign(si3[k-1])
    # ss4[k] = si4[k] + alpha1*si4[k-1] + alpha2*np.sign(si4[k-1])

    if k == 1:
        ss1[1] =0
        ss2[1] =0
        ss3[1] =0
        ss4[1] =0

    else:
        ss1[k+1] = si1[k+1] + alpha1 * si1[k]
        ss2[k+1] = si2[k+1] + alpha1 * si2[k]
        ss3[k+1] = si3[k+1] + alpha1 * si3[k]
        ss4[k+1] = si4[k+1] + alpha1 * si4[k]



    if k == 1:
            u1[1] = 0
            u2[1] = 0
            u3[1] = 0
            u4[1] = 0
    else:
        mfa1[k] = mfa1[k - 1] + (rho * phi1[k]) / (lambda_ + abs(phi1[k])**2) * si1[k]
        mfa2[k] = mfa2[k - 1] + (rho * phi2[k]) / (lambda_ + abs(phi2[k])**2) * si2[k]
        mfa3[k] = mfa3[k - 1] + (rho * phi3[k]) / (lambda_ + abs(phi3[k])**2) * si3[k]
        mfa4[k] = mfa4[k - 1] + (rho * phi4[k]) / (lambda_ + abs(phi4[k])**2) * si4[k]


    if k == 1:
        sm1[1] = 0
        sm2[1] = 0
        sm3[1] = 0
        sm4[1] = 0
    else:
        sm1[k] = sm1[k-1] + (ss1[k]-alpha1*si1[k]-epsilon*T*np.sign(ss1[k])-si1[k]) / (phi1[k])
        sm1[k] = sm1[k-1] + (ss1[k]-alpha1*si1[k]-epsilon*T*np.sign(ss1[k])-si1[k]) / (phi1[k])
        sm1[k] = sm1[k-1] + (ss1[k]-alpha1*si1[k]-epsilon*T*np.sign(ss1[k])-si1[k]) / (phi1[k])
        sm1[k] = sm1[k-1] + (ss1[k]-alpha1*si1[k]-epsilon*T*np.sign(ss1[k])-si1[k]) / (phi1[k])

        # sm1[k] = sm1[k-1] + (-beta * si1[k] - epsilon * T * np.sign(ss1[k])-alpha1 * si1[k] - si1[k]) / (phi1[k])
        # sm2[k] = sm2[k-1] + (-beta * si2[k] - epsilon * T * np.sign(ss2[k])-alpha1 * si2[k] - si2[k]) / (phi2[k])
        # sm3[k] = sm3[k-1] + (-beta * si3[k] - epsilon * T * np.sign(ss3[k])-alpha1 * si3[k] - si3[k]) / (phi3[k])
        # sm4[k] = sm4[k-1] + (-beta * si4[k] - epsilon * T * np.sign(ss4[k])-alpha1 * si4[k] - si4[k]) / (phi4[k])

    
    if k == 1:
        u1[1] = 0.1
        u2[1] = 0.1
        u3[1] = 0.1
        u4[1] = 0.1
    else:
        u1[k] = mfa1[k] + gamma1 * sm1[k]
        u2[k] = mfa2[k] + gamma2 * sm2[k]
        u3[k] = mfa3[k] + gamma3 * sm3[k]
        u4[k] = mfa4[k] + gamma4 * sm4[k]



    y1[1] = 0.1
    y2[1] = 0.1
    y3[1] = 0.1
    y4[1] = 0.1


    y1[k + 1] = m / (rT * 0.1) * u1[k]
    y2[k + 1] = m / (rT * 0.1) * u2[k]
    y3[k + 1] = m / (rT * 0.3) * u3[k]
    y4[k + 1] = m / (rT * 0.3) * u4[k]

    e1[k] = abs(y1[k]-w5[k])
    e2[k] = abs(y2[k]-y1[k])
    e3[k] = abs(y3[k]-y2[k])
    e4[k] = abs(y4[k]-w6[k])



# plt.figure(figsize=(8, 6))


def plot_output_figure():
    plt.figure()
    plt.plot(y1[:-1], '-k', markersize=4, label='Y1')
    plt.plot(y2[:-1], '-', markersize=4, label='Y2')
    plt.plot(y3[:-1], '--g', markersize=4, label='Y3')
    plt.plot(y4[:-1], '--b', markersize=4, label='Y4')
    # plt.plot(yd[:-1], '-*g', markersize=4, label='Y2')
    plt.plot(w5[:-1], '-*r', markersize=2, label='w5')
    plt.plot(w6[:-1], '-*r', markersize=2, label='w5')
    plt.grid(True)
    plt.xlabel("time step", fontsize=12)
    plt.ylabel("outputs", fontsize=12)
    plt.title("Without Second Order Sliding Dynamics",fontsize=10)
    plt.legend(fontsize=10)




# Plotting containment errors
# plt.figure(figsize=(8, 6))
# plt.plot(e1[:-1], '-r', markersize=4, label='Error 1 (E1)')
# plt.plot(e2[:-1], '-b', markersize=4, label='Error 2 (E2)')
# plt.plot(e3[:-1], '--g', markersize=4, label='Error 3 (E3)')
# plt.plot(e4[:-1], '--k', markersize=4, label='Error 4 (E4)')
# plt.title('Containment Errors')
# plt.xlabel('Time Step')
# plt.ylabel('Error Magnitude')
# plt.legend()
# plt.grid()
# plt.show()

def plot_figure():

    plt.figure()
    plt.plot(phi1[:-1], '-k', markersize=4, label='phi1')
    plt.plot(phi2[:-1], '-', markersize=4, label='phi2')
    plt.plot(phi3[:-1], '--g', markersize=4, label='phi3')
    plt.plot(phi4[:-1], '--b', markersize=4, label='phi4')
    # plt.plot(yd[:-1], '-*g', markersize=4, label='Y2')
    # plt.plot(w5[:-1], '-*r', markersize=2, label='w5')
    # plt.plot(w6[:-1], '-*r', markersize=2, label='w5')
    plt.grid(True)
    plt.xlabel("time step", fontsize=12)
    plt.ylabel("outputs", fontsize=12)
    plt.title("Without Second Order Sliding Dynamics",fontsize=10)





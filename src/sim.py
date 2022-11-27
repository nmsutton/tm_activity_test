# Create time points
Tmin, Tmax, dt = 0, 1000, 0.025  # Step size
T = np.arange(Tmin, Tmax + dt, dt)

# Declare the current to be used
I_ext = 150

## Nine parameters for IM model
# Axo_Axonic
k, a, b, d, C, vr, vt, c, vpeak = 3.961462878, 0.004638608, 8.683644937, 15, 165, -57.09978287, \
                                  -51.71875628, -73.96850421, 27.79863559

preNeuronType = "Axo-Axonic"

## Make a state vector that has a (v, u) pair for each timestep
s = np.zeros((len(T), 2))

# Create a vector to store spike times
spike_times = np.array([])
isi_mode_list = np.array([])


## Initial values
s[0, 0] = vr
s[0, 1] = 0


# Note that s1[0] is v, s1[1] is u. This is Izhikevich equation in vector form
def s_dt(s1, I):
  v_dt = (k*(s1[0] - vr)*(s1[0] - vt) - s1[1] + I)*(1/C)
  u_dt = a*(b*(s1[0] - vr) - s1[1])
  return np.array([v_dt, u_dt])


## SIMULATE
for t in range(len(T)-1):
  # Calculate the four constants of Runge-Kutta method
  k_1 = s_dt(s[t], I_ext)
  k_2 = s_dt(s[t] + 0.5*dt*k_1, I_ext)
  k_3 = s_dt(s[t] + 0.5*dt*k_2, I_ext)
  k_4 = s_dt(s[t] + dt*k_3, I_ext)

  s[t+1] = s[t] + (1.0/6)*dt*(k_1 + 2*k_2 + 2*k_3 + k_4)

  # Reset the neuron if it has spiked
  if s[t+1, 0] >= vpeak:
    s[t, 0]   = vpeak # add Dirac pulse for visualisation
    s[t+1, 0] = c   # Reset to resting potential
    s[t+1, 1] += d  # Update recovery variable
    spike_times = np.append(spike_times, math.ceil(t*dt))

v = s[:, 0]
u = s[:, 1]


## Nine parameters for IM model
# Pyramidal
postNeuronType = "Pyramidal"
a, b, c, d, k, vr, vt, C, vpeak = 0.00838350334098279, -42.5524776883928, -38.8680990294091, 588.0, \
                                  0.792338703789581, -63.2044008171655, -33.6041733124267, 366.0, 35.8614648558726

# Synaptic Parameters and equations
synaptic_event_times = list(spike_times)

# TPM parameters
g0, tau_d, tau_f, tau_r, Utilization, e_rev = 3.644261648, 10.71107251, 17.20004939, 435.8103009, 0.259914361, -70

g0 /= Utilization
tm_model = 'Keivan'


def synaptic_event(delta_t_, g0, tau_d_, tau_r_, tau_f_, utilization, x0, y0, u0):
    if tm_model != 'Carlsim':
        # TM Model that depends on tau_d
        tau1r = tau_d_ / ((tau_d_ - tau_r_) if tau_d_ != tau_r_ else 1e-13)
        y_ = y0 * math.exp(-delta_t_ / tau_d_)
        x_ = 1 + (x0 - 1 + tau1r * y0) * math.exp(-delta_t_ / tau_r_) - tau1r * y_
        u_ = u0 * math.exp(-delta_t_ / tau_f_)
        u0 = u_ + utilization * (1 - u_)
        y0 = y_ + u0 * x_
        x0 = x_ - u0 * x_
        g = g0 * y0
    else:
        # Carlsim's TM Model
        x_ = 1 + (x0 - 1) * math.exp(-delta_t_ / tau_r_)
        u_ = utilization + (u0 - utilization) * math.exp(-delta_t_ / tau_f_)
#         print(f'x0={x0:.3f}, x_={x_:.3f}, u0={u0:.3f}, u_={u_:.3f}')
        g = g0 * u_ * x_
        x0 = x_ - u_ * x_
        u0 = u_ + utilization * (1 - u_)
    return g, x0, y0, u0


def synaptic_current(g, delta_t_, tau_d_, e_syn_):
    return g * math.exp(-delta_t_ / tau_d_) * e_syn_


# Initialize synaptic state variables
conductance, x0, y0 = 0.0, 1.0, 0.0
I = 0.0
if tm_model == 'Carlsim':
    u0 = Utilization
else:
    u0 = 0

# Make a state vector that has a (v, u) pair for each timestep
s = np.zeros((len(T), 2))

# Initial Izhikevich state variables
s[0, 0] = vr
s[0, 1] = 0


# Note that s1[0] is v, s1[1] is u. This is Izhikevich equation in vector form
def s_dt(s1, I):
    v_dt = (k * (s1[0] - vr) * (s1[0] - vt) - s1[1] + I) * (1 / C)
    u_dt = a * (b * (s1[0] - vr) - s1[1])
    return np.array([v_dt, u_dt])

# SIMULATE
next_synaptic_event_time, delta_t, I_syn = synaptic_event_times[0], 0.0, [0]
synaptic_event_time = next_synaptic_event_time
spike_times = np.array([])
for t in range(len(T) - 1):
    v = s[t, 0]
    e_syn = v - e_rev
    time = T[t]

    if next_synaptic_event_time <= time:
        inter_event_time = next_synaptic_event_time - synaptic_event_time
        conductance, x0, y0, u0 = synaptic_event(inter_event_time, g0, tau_d, tau_r, tau_f, Utilization, x0, y0, u0)
        # print(conductance, time, synaptic_event_time, next_synaptic_event_time)
        synaptic_event_time = time
        if len(synaptic_event_times) > 1:
            del synaptic_event_times[0]
            next_synaptic_event_time = synaptic_event_times[0]
        else:
            next_synaptic_event_time = math.inf

    delta_t = time - synaptic_event_time
    if delta_t >= 0:
        I = synaptic_current(conductance, delta_t, tau_d, e_syn)
    I_syn.append(I)

    # Calculate the four constants of Runge-Kutta method
    k_1 = s_dt(s[t], -I)
    k_2 = s_dt(s[t] + 0.5 * dt * k_1, -I)
    k_3 = s_dt(s[t] + 0.5 * dt * k_2, -I)
    k_4 = s_dt(s[t] + dt * k_3, -I)

    s[t + 1] = s[t] + (1.0 / 6) * dt * (k_1 + 2 * k_2 + 2 * k_3 + k_4)

    # Reset the neuron if it has spiked
    if s[t + 1, 0] >= vpeak:
        s[t, 0] = vpeak  # add Dirac pulse for visualisation
        s[t + 1, 0] = c  # Reset to resting potential
        s[t + 1, 1] += d  # Update recovery variable
        spike_times = np.append(spike_times, math.floor(t * dt))

v = s[:, 0]
u = s[:, 1]

# First define function to flip the sign of the current
def sign(lst): 
    return [ -i for i in lst ]

# Compare to CARLsim simulation output
FH = np.loadtxt("HC_IM_05_26_aac_pyr_I_150pA_fast_1_slow_0.txt")
I = FH[1::2]
V = FH[0::2]
I = sign(I)
I = I[0:len(I_syn)-1]
V = V[0:len(I_syn)-1]
ax1 = plt.subplot(211)
ax1.plot(T,np.append(V, [V[-1]]), label = "CARLsim TM") # added the last value of V to ensure the same length as time vector
plt.ylabel('Membrane potential (mV)')
ax2 = plt.subplot(212)
ax2.plot(T,np.append(I, [I[-1]]), label = "CARLsim TM") # added the last value of I to ensure the same length as time vector
plt.ylabel('Synaptic Current (pA)')
plt.xlabel('Time (ms)')
plt.legend(loc = "center right")
plt.tight_layout()

## Plot the membrane potential
ax1 = plt.subplot(211)
ax1.plot(T, v, color = "orange", linestyle='dotted', label = "Keivan TM", alpha=0.85)
plt.ylabel('Membrane potential (mV)')
plt.title(f"{postNeuronType}")
ax2 = plt.subplot(212)
ax2.plot(T, I_syn, color = "orange", linestyle='dotted', label = "Keivan TM", alpha=0.85)
plt.ylabel('Synaptic Current (pA)')
plt.xlabel('Time (ms)')
plt.legend(loc = "center right")
plt.tight_layout()
fileOutputName = preNeuronType + '_' + postNeuronType + '_' + str(I_ext) + 'pA' + \
      '_CARLsim_vs_Keivan_superimposed.png'
plt.savefig(fileOutputName, dpi=800)
plt.clf()


# Look at the error between CARLsim and python computed synaptic signal
I = np.append(I, [I[-1]])
V = np.append(V, [V[-1]])

af = scipy.fft.fft(I_syn)
bf = scipy.fft.fft(I)
c = scipy.ifft(af * scipy.conj(bf))
time_shift = np.argmax(abs(c))
#         print(time_shift)

I_2 = I[:-1]
I_2 = I_2[::int(1/dt)]
I_syn_2 = I_syn[0+time_shift:]
I_syn_2 = I_syn_2[::int(1/dt)]
I_syn_2 = np.array(np.array(I_syn_2,dtype=np.float32))
V_2 = V[:-1]
V_2 = V_2[::int(1/dt)]
v_2 = v[0+time_shift:]
v_2 = v_2[::int(1/dt)]

if len(V_2) == len(v_2):
    ax1 = plt.subplot(211)
    ax1.plot(abs(V_2 - v_2))
    plt.ylabel('Membrane potential (mV)')
    plt.title(f"{postNeuronType}")
    ax2 = plt.subplot(212)
    ax2.plot(abs(I_2 - I_syn_2))
    plt.ylabel('I_syn diff (HCO - Carlsim) (pA)')
    plt.xlabel('Time (ms)')
    plt.tight_layout()
    fileOutputName = preNeuronType + '_' + postNeuronType + '_' + str(I_ext) + 'pA' + \
          '_CARLsim_vs_Keivan_error.png'
    plt.savefig(fileOutputName, dpi=800)
    plt.clf()

    # Append min, max, mean, and median of errors between V and I
    pctErrorV = abs((V_2 - v_2)/v_2)
    pctErrorI = abs((I_2 - I_syn_2)/I_syn_2)
    maxErrorV = max(pctErrorV[~np.isnan(pctErrorV)])
    maxErrorI = max(pctErrorI[~np.isnan(pctErrorI)])
    minErrorV = min(pctErrorV[~np.isnan(pctErrorV)])
    minErrorI = min(pctErrorI[~np.isnan(pctErrorI)])
    meanErrorV = np.mean(pctErrorV[~np.isnan(pctErrorV)])
    meanErrorI = np.mean(pctErrorI[~np.isnan(pctErrorI)])
    medianErrorV = np.median(pctErrorV[~np.isnan(pctErrorV)])
    medianErrorI = np.median(pctErrorI[~np.isnan(pctErrorI)])
    mismatchV = sum(abs(V_2-v_2)*dt)/sum(abs((V_2 + v_2)/2))
    mismatchI = sum(abs(I_2-I_syn_2)*dt)/sum(abs((I_2 + I_syn_2)/2))

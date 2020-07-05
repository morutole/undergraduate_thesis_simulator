import numpy as np
from numpy.core.defchararray import index
import pandas as pd
import matplotlib.pyplot as plt

def read_raw_csv():
    raw_simulate_df = pd.read_csv("sample.csv")
    search_word = "sat_position"
    raw_position_df = raw_simulate_df.loc[:,raw_simulate_df.columns.str.contains(search_word)]
    search_word = "sat_velocity_i" 
    raw_velocity_df = raw_simulate_df.loc[:,raw_simulate_df.columns.str.contains(search_word)]

    return raw_position_df.values, raw_velocity_df.values

def read_M_csv():
    raw_M_df = pd.read_csv("ConvMatrix.csv", header = None)
    raw_M_nd = raw_M_df.values
    return raw_M_nd

def read_general_csv(file_name):
    value_df = pd.read_csv(file_name, header = None)
    value_nd = value_df.values
    position_nd = value_nd[:, 0:3]
    velocity_nd = value_nd[:, 3:6]

    return position_nd, velocity_nd

def extract_observed_correspondent_value(GPS_period, true_position_nd, true_velocity_nd):
    true_time_scale = np.linspace(0, simulation_time, true_position_nd.shape[0]).astype(int)
    ind = np.where(true_time_scale%GPS_period == 0)[0]
    true_observed_position_nd = true_position_nd[ind]
    true_observed_velocity_nd = true_velocity_nd[ind]

    return true_observed_position_nd, true_observed_velocity_nd

def get_label(index, kind):
    if index == 0:
        res = "x"
    elif index == 1:
        res = "y"
    elif index == 2:
        res = "z"
    
    if kind == "velocity":
        return "V" + res + "[m/s]"
    else:
        return res + "[m]"

raw_position_nd, raw_velocity_nd = read_raw_csv()
true_position_nd, true_velocity_nd = read_general_csv("true.csv")
observed_position_nd, observed_velocity_nd = read_general_csv("observed.csv")
estimate_position_nd, estimate_velocity_nd = read_general_csv("estimate.csv")
Conv_M_nd = read_M_csv()

simulation_time = true_position_nd.shape[0] - 1 #開始時点まで含んでいるため
raw_time_scale = np.linspace(0, simulation_time, raw_position_nd.shape[0])
assert simulation_time%(observed_position_nd.shape[0] - 1) == 0, "GPS_periodがおかしい"
GPS_period = simulation_time/(observed_position_nd.shape[0] - 1)

true_observed_position_nd, true_observed_velocity_nd = extract_observed_correspondent_value(GPS_period, true_position_nd, true_velocity_nd)

observed_time_scale = np.linspace(0, simulation_time, observed_position_nd.shape[0])

def plot_absolute_graph():
    for i in range(3):
        plt.plot(raw_time_scale, raw_position_nd[:, i], label = "true")
        plt.plot(raw_time_scale, estimate_position_nd[:, i], label = "estimate")
        plt.xlabel("time[s]")
        plt.ylabel(get_label(i, "position"))
        plt.legend(loc = "upper right")
        plt.savefig("abs_x" + str(i) + ".png",bbox_inches = "tight", pad_inches = 0.05)
        plt.show()

    for i in range(3):
        plt.plot(raw_time_scale, raw_velocity_nd[:, i], label = "true")
        plt.plot(raw_time_scale, estimate_velocity_nd[:, i], label = "estimate")
        plt.xlabel("time[s]")
        plt.ylabel(get_label(i, "velocity"))
        plt.legend(loc = "upper right")
        plt.savefig("abs_v" + str(i) + ".png",bbox_inches = "tight", pad_inches = 0.05)
        plt.show()

def plot_difference_graph():
    for i in range(3):
        plt.scatter(observed_time_scale, observed_position_nd[:, i] - true_observed_position_nd[:, i], label = "observed", color = "C0")
        plt.scatter(raw_time_scale, estimate_position_nd[:, i] - raw_position_nd[:, i], label = "estimate", color = "C1")
        plt.xlabel("time[s]")
        plt.ylabel(get_label(i, "position"))
        plt.ylim(-2.0, 2.0)
        plt.legend(loc = "upper right")
        plt.savefig("diff_x" + str(i) + ".png",bbox_inches = "tight", pad_inches = 0.05)
        plt.show()

    for i in range(3):
        plt.scatter(observed_time_scale, observed_velocity_nd[:, i] - true_observed_velocity_nd[:, i], label = "observed", color = "C0")
        plt.scatter(raw_time_scale, estimate_velocity_nd[:, i] - raw_velocity_nd[:, i], label = "estimate", color = "C1")
        plt.xlabel("time[s]")
        plt.ylabel(get_label(i, "velocity"))
        plt.ylim(-0.2, 0.2)
        plt.legend(loc = "upper right")
        plt.savefig("diff_v" + str(i) + ".png",bbox_inches = "tight", pad_inches = 0.05)
        plt.show()

def plot_rms_graph():
    plt.scatter(observed_time_scale, np.linalg.norm(observed_position_nd - true_observed_position_nd, axis = 1), label = "observed", color = "C0")
    plt.scatter(raw_time_scale, np.linalg.norm(estimate_position_nd - raw_position_nd, axis = 1), label = "estimate", color = "C1")
    plt.xlabel("time[s]")
    plt.ylabel("position[m]")
    plt.ylim(0, 2.5)
    plt.legend(loc = "upper right")
    plt.savefig("rms_r.png", bbox_inches = "tight", pad_inches = 0.05)
    plt.show()

    plt.scatter(observed_time_scale, np.linalg.norm(observed_velocity_nd - true_observed_velocity_nd, axis = 1), label = "observed", color = "C0")
    plt.scatter(raw_time_scale, np.linalg.norm(estimate_velocity_nd - raw_velocity_nd, axis = 1), label = "estimate", color = "C1")
    plt.xlabel("time[s]")
    plt.ylabel("velocity[m/s]")
    plt.ylim(0, 0.25)
    plt.legend(loc = "upper right")
    plt.savefig("rms_v.png", bbox_inches = "tight", pad_inches = 0.05)
    plt.show()

def plot_M_graph(t_max):
    index_sup = t_max*1000+1
    time_scale = np.linspace(0,t_max,index_sup) #1e-3秒
    
    plt.plot(time_scale, np.sqrt(Conv_M_nd[0:index_sup, 0]))
    plt.xlabel("time[s]")
    plt.ylabel(r"$\sigma_x[\mathrm{m}]$")
    plt.ylim(0.0, 5.5)
    plt.savefig("sigma_x.png", bbox_inches = "tight", pad_inches = 0.05)
    plt.show()
    plt.plot(time_scale, np.sqrt(Conv_M_nd[0:index_sup, 1]))
    plt.xlabel("time[s]")
    plt.ylabel(r"$\sigma_x[\mathrm{m}]$")
    plt.ylim(0.0, 5.5)
    plt.savefig("sigma_y.png", bbox_inches = "tight", pad_inches = 0.05)
    plt.show()
    plt.plot(time_scale, np.sqrt(Conv_M_nd[0:index_sup, 2]))
    plt.xlabel("time[s]")
    plt.ylabel(r"$\sigma_x[\mathrm{m}]$")
    plt.ylim(0.0, 5.5)
    plt.savefig("sigma_z.png", bbox_inches = "tight", pad_inches = 0.05)
    plt.show()

    plt.plot(time_scale, np.sqrt(Conv_M_nd[0:index_sup, 3]))
    plt.xlabel("time[s]")
    plt.ylabel(r"$\sigma_{vx}[\mathrm{m/s}]$")
    plt.ylim(0.0, 0.55)
    plt.savefig("sigma_vx.png", bbox_inches = "tight", pad_inches = 0.05)
    plt.show()
    plt.plot(time_scale, np.sqrt(Conv_M_nd[0:index_sup, 4]))
    plt.xlabel("time[s]")
    plt.ylabel(r"$\sigma_{vy}[\mathrm{m/s}]$")
    plt.ylim(0.0, 0.55)
    plt.savefig("sigma_vy.png", bbox_inches = "tight", pad_inches = 0.05)
    plt.show()
    plt.plot(time_scale, np.sqrt(Conv_M_nd[0:index_sup, 5]))
    plt.xlabel("time[s]")
    plt.ylabel(r"$\sigma_{vz}[\mathrm{m/s}]$")
    plt.ylim(0.0, 0.55)
    plt.savefig("sigma_vz.png", bbox_inches = "tight", pad_inches = 0.05)
    plt.show()

    plt.plot(time_scale, np.sqrt(Conv_M_nd[0:index_sup, 6]))
    plt.xlabel("time[s]")
    plt.ylabel(r"$\sigma_{ax}[\mathrm{m/s^2}]$")
    plt.ylim(0.0, 0.11)
    plt.savefig("sigma_ax.png", bbox_inches = "tight", pad_inches = 0.05)
    plt.show()
    plt.plot(time_scale, np.sqrt(Conv_M_nd[0:index_sup, 7]))
    plt.xlabel("time[s]")
    plt.ylabel(r"$\sigma_{ay}[\mathrm{m/s^2}]$")
    plt.ylim(0.0, 0.11)
    plt.savefig("sigma_ay.png", bbox_inches = "tight", pad_inches = 0.05)
    plt.show()
    plt.plot(time_scale, np.sqrt(Conv_M_nd[0:index_sup, 8]))
    plt.xlabel("time[s]")
    plt.ylabel(r"$\sigma_{az}[\mathrm{m/s^2}]$")
    plt.ylim(0.0, 0.11)
    plt.savefig("sigma_az.png", bbox_inches = "tight", pad_inches = 0.05)
    plt.show()

plot_absolute_graph()
plot_difference_graph()
plot_rms_graph()
plot_M_graph(6)

print(np.mean(np.linalg.norm(estimate_position_nd - raw_position_nd, axis = 1)))
print(np.mean(np.linalg.norm(observed_position_nd - true_observed_position_nd, axis = 1)))
print(np.mean(np.linalg.norm(estimate_velocity_nd - raw_velocity_nd, axis = 1)))
print(np.mean(np.linalg.norm(observed_velocity_nd - true_observed_velocity_nd, axis = 1)))